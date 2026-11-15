//---------------------------------------------------------------------------
// SXMesh.cpp — VC++ · C++11 판. 머리말은 SXMesh.h 참고.
//
// NMF 형식(텍스트, INI 비슷함)
//   [Info]                     FormatVersion, nNode, nElement, nTexturePoint, iNodeShape, GLAlpha,
//                              nMTLData, PropertyBit, PointSize, LineWidth,
//                              SXNNode=이름/행/열,...   SXNElement=이름/행/열,...
//   [pMTLData(m)]              m = 0..nMTLData-1. ID, Name, map_Kd, map_Ks, map_d, bump, Ka, Kd, Ks,
//                              Shininess, GLAlpha, 그 밖의 키는 그대로
//   [그 밖의 섹션]             예: [ManikinOBB]. 키=값 을 그대로 보존
//   [변수 이름]                행마다 한 줄, 값은 쉼표로 구분. 행의 값 개수 = 그 줄의 값 개수
//
// 옛 C++Builder 판 파일과의 차이(읽기는 둘 다 된다)
//   - 표시 옵션 키: 옛 판은 bData=<PropertyBit 정수> 로 썼다. 새 판은 PropertyBit=. 읽을 때 둘 다 본다.
//   - 값 정밀도: 옛 판은 유효숫자 3자리(FloatToStrF). 새 판은 왕복 정밀도.
//   - 1열 변수: 옛 판은 한 줄에 모두 쓰고 끝에 쉼표. 새 판은 다른 변수처럼 한 줄에 한 행.
//   - 재질 섹션 제목: 옛 판은 pMTLData(ID). 새 판은 pMTLData(순번) + ID= 키. 읽을 때 둘 다 된다.
//   - 문자열: 옛 판은 CP949. 새 판은 UTF-8. 읽을 때 줄마다 UTF-8 이 아니면 CP949 로 읽는다(Windows).
//---------------------------------------------------------------------------
#include "SXMesh.h"

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cwchar>
#include <exception>
#include <fstream>
#include <istream>
#include <iterator>
#include <map>
#include <new>
#include <ostream>
#include <set>
#include <stdexcept>

#ifdef _WIN32
	#ifndef NOMINMAX
		#define NOMINMAX
	#endif
	#ifndef WIN32_LEAN_AND_MEAN
		#define WIN32_LEAN_AND_MEAN
	#endif
	#include <windows.h>
#endif

#ifdef _MSC_VER
	#include <locale.h>
#endif

namespace {

//===========================================================================
// 문자 부호화
//===========================================================================
void AppendCodePoint(std::wstring& out, unsigned long cp)
{
#if WCHAR_MAX <= 0xFFFF
	if (cp >= 0x10000UL) {   // UTF-16(Windows): 대리 쌍
		cp -= 0x10000UL;
		out.push_back(static_cast<wchar_t>(0xD800UL + (cp >> 10)));
		out.push_back(static_cast<wchar_t>(0xDC00UL + (cp & 0x3FFUL)));
		return;
	}
#endif
	out.push_back(static_cast<wchar_t>(cp));
}

bool DecodeUtf8(const char* s, std::size_t n, std::wstring& out)
{
	out.clear();
	out.reserve(n);
	std::size_t i = 0;
	while (i < n) {
		const unsigned char c = static_cast<unsigned char>(s[i]);
		unsigned long cp;
		int extra;
		if (c < 0x80)                { cp = c;        extra = 0; }
		else if ((c & 0xE0) == 0xC0) { cp = c & 0x1F; extra = 1; }
		else if ((c & 0xF0) == 0xE0) { cp = c & 0x0F; extra = 2; }
		else if ((c & 0xF8) == 0xF0) { cp = c & 0x07; extra = 3; }
		else return false;
		if (extra > 0 && i + static_cast<std::size_t>(extra) >= n) return false;   // 뒤 바이트가 모자람
		for (int k = 1; k <= extra; ++k) {
			const unsigned char cc = static_cast<unsigned char>(s[i + k]);
			if ((cc & 0xC0) != 0x80) return false;
			cp = (cp << 6) | (cc & 0x3F);
		}
		if ((extra == 1 && cp < 0x80) || (extra == 2 && cp < 0x800) || (extra == 3 && cp < 0x10000) ||
		    cp > 0x10FFFFUL || (cp >= 0xD800 && cp <= 0xDFFF))
			return false;
		AppendCodePoint(out, cp);
		i += static_cast<std::size_t>(extra) + 1;
	}
	return true;
}

std::string EncodeUtf8(const std::wstring& w)
{
	std::string out;
	out.reserve(w.size());
	for (std::size_t i = 0; i < w.size(); ++i) {
		unsigned long cp = static_cast<unsigned long>(w[i]);
#if WCHAR_MAX <= 0xFFFF
		cp &= 0xFFFFUL;
		if (cp >= 0xD800 && cp <= 0xDBFF && i + 1 < w.size()) {   // 대리 쌍
			const unsigned long lo = static_cast<unsigned long>(w[i + 1]) & 0xFFFFUL;
			if (lo >= 0xDC00 && lo <= 0xDFFF) {
				cp = 0x10000UL + ((cp - 0xD800UL) << 10) + (lo - 0xDC00UL);
				++i;
			}
		}
#endif
		if (cp >= 0xD800 && cp <= 0xDFFF) cp = 0xFFFD;   // 짝 없는 대리 문자
		if (cp > 0x10FFFFUL) cp = 0xFFFD;
		if (cp < 0x80) {
			out.push_back(static_cast<char>(cp));
		} else if (cp < 0x800) {
			out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
			out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
		} else if (cp < 0x10000) {
			out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
			out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
			out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
		} else {
			out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
			out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
			out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
			out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
		}
	}
	return out;
}

// UTF-8 이면 UTF-8 로, 아니면 옛 파일(CP949)로 읽는다.
std::wstring DecodeText(const char* s, std::size_t n)
{
	std::wstring w;
	if (DecodeUtf8(s, n, w)) return w;
#ifdef _WIN32
	if (n > 0 && n <= static_cast<std::size_t>(INT_MAX)) {
		const UINT pages[2] = { 949, CP_ACP };
		for (int k = 0; k < 2; ++k) {
			const DWORD flags = (k == 0) ? MB_ERR_INVALID_CHARS : 0;
			const int len = MultiByteToWideChar(pages[k], flags, s, static_cast<int>(n), 0, 0);
			if (len > 0) {
				w.assign(static_cast<std::size_t>(len), L'\0');
				MultiByteToWideChar(pages[k], flags, s, static_cast<int>(n), &w[0], len);
				return w;
			}
		}
	}
#endif
	w.clear();   // 마지막 수단: Latin-1
	for (std::size_t i = 0; i < n; ++i) w.push_back(static_cast<wchar_t>(static_cast<unsigned char>(s[i])));
	return w;
}

//===========================================================================
// 숫자 (프로그램이 setlocale 을 바꿔도 소수점은 늘 '.')
//===========================================================================
#ifdef _MSC_VER
_locale_t CLocale()
{
	static _locale_t loc = _create_locale(LC_NUMERIC, "C");
	return loc;
}
#endif

bool ParseDouble(const char* b, const char* e, double& v)
{
	const std::size_t len = static_cast<std::size_t>(e - b);
	if (len == 0) return false;
	char small[128];
	std::string big;
	const char* p;
	if (len < sizeof(small)) { std::memcpy(small, b, len); small[len] = '\0'; p = small; }
	else { big.assign(b, e); p = big.c_str(); }
	char* endp = 0;
#ifdef _MSC_VER
	v = _strtod_l(p, &endp, CLocale());
#else
	v = std::strtod(p, &endp);
#endif
	return endp == p + len;
}

bool ParseFloat(const char* b, const char* e, float& f)
{
	double d;
	if (!ParseDouble(b, e, d)) return false;
	f = static_cast<float>(d);
	return true;
}

bool ParseLong(const char* b, const char* e, long long& v)
{
	const std::size_t len = static_cast<std::size_t>(e - b);
	if (len == 0 || len >= 64) return false;
	char buf[64];
	std::memcpy(buf, b, len);
	buf[len] = '\0';
	char* endp = 0;
	errno = 0;
	v = std::strtoll(buf, &endp, 10);
	if (endp == buf + len && errno == 0) return true;
	double d;   // 옛 파일이 정수를 "3.000" 처럼 썼을 경우
	if (ParseDouble(b, e, d) && d >= -9.0e18 && d <= 9.0e18 && d == std::floor(d)) {
		v = static_cast<long long>(d);
		return true;
	}
	return false;
}

bool ParseInt(const char* b, const char* e, int& v)
{
	long long t;
	if (!ParseLong(b, e, t) || t < INT_MIN || t > INT_MAX) return false;
	v = static_cast<int>(t);
	return true;
}

void AppendFormatted(std::string& out, const char* fmt, double v)
{
	char buf[64];
#ifdef _MSC_VER
	_sprintf_s_l(buf, sizeof(buf), fmt, CLocale(), v);
#else
	std::snprintf(buf, sizeof(buf), fmt, v);
#endif
	out += buf;
}

// 다시 읽었을 때 같은 값이 되는 가장 짧은 쪽(유효숫자 15 → 17)
void AppendDouble(std::string& out, double v)
{
	char buf[64];
#ifdef _MSC_VER
	_sprintf_s_l(buf, sizeof(buf), "%.15g", CLocale(), v);
#else
	std::snprintf(buf, sizeof(buf), "%.15g", v);
#endif
	double back;
	if (ParseDouble(buf, buf + std::strlen(buf), back) && back == v) { out += buf; return; }
	AppendFormatted(out, "%.17g", v);
}

void AppendFloat(std::string& out, float f)
{
	char buf[64];
#ifdef _MSC_VER
	_sprintf_s_l(buf, sizeof(buf), "%.7g", CLocale(), static_cast<double>(f));
#else
	std::snprintf(buf, sizeof(buf), "%.7g", static_cast<double>(f));
#endif
	float back;
	if (ParseFloat(buf, buf + std::strlen(buf), back) && back == f) { out += buf; return; }
	AppendFormatted(out, "%.9g", static_cast<double>(f));
}

void AppendLong(std::string& out, long long v)
{
	char buf[32];
	std::snprintf(buf, sizeof(buf), "%lld", v);
	out += buf;
}

void AppendValue(std::string& out, double v) { AppendDouble(out, v); }
void AppendValue(std::string& out, int v)    { AppendLong(out, v); }

bool ParseValue(const char* b, const char* e, double& v) { return ParseDouble(b, e, v); }
bool ParseValue(const char* b, const char* e, int& v)    { return ParseInt(b, e, v); }

//===========================================================================
// 줄 · 토큰
//===========================================================================
struct Line {
	const char* b;
	const char* e;
};

struct Token {
	const char* b;
	const char* e;
};

bool IsSpace(char c) { return c == ' ' || c == '\t'; }

void SplitLines(const char* buf, std::size_t size, std::vector<Line>& lines)
{
	const char* p = buf;
	const char* end = buf + size;
	if (size >= 3 && static_cast<unsigned char>(p[0]) == 0xEF && static_cast<unsigned char>(p[1]) == 0xBB &&
	    static_cast<unsigned char>(p[2]) == 0xBF)
		p += 3;   // UTF-8 BOM
	while (p < end) {
		const char* nl = static_cast<const char*>(std::memchr(p, '\n', static_cast<std::size_t>(end - p)));
		const char* le = nl ? nl : end;
		Line ln;
		ln.b = p;
		ln.e = le;
		while (ln.e > ln.b && (ln.e[-1] == '\r')) --ln.e;
		lines.push_back(ln);
		p = nl ? nl + 1 : end;
	}
}

// 쉼표로 나눈다. 끝의 빈 토큰(옛 판의 끝 쉼표)은 버린다. 가운데 빈 토큰이 있으면 false.
bool Tokenize(const Line& ln, std::vector<Token>& toks)
{
	toks.clear();
	const char* p = ln.b;
	while (p < ln.e && IsSpace(*p)) ++p;
	if (p == ln.e) return true;   // 빈 줄 = 값 0개
	const char* s = ln.b;
	for (const char* q = ln.b; ; ++q) {
		if (q == ln.e || *q == ',') {
			Token t;
			t.b = s;
			t.e = q;
			while (t.b < t.e && IsSpace(*t.b)) ++t.b;
			while (t.e > t.b && IsSpace(t.e[-1])) --t.e;
			toks.push_back(t);
			if (q == ln.e) break;
			s = q + 1;
		}
	}
	while (!toks.empty() && toks.back().b == toks.back().e) toks.pop_back();
	for (std::size_t i = 0; i < toks.size(); ++i)
		if (toks[i].b == toks[i].e) return false;
	return true;
}

std::wstring Trim(const std::wstring& s)
{
	std::size_t b = 0, e = s.size();
	while (b < e && (s[b] == L' ' || s[b] == L'\t')) ++b;
	while (e > b && (s[e - 1] == L' ' || s[e - 1] == L'\t')) --e;
	return s.substr(b, e - b);
}

std::vector<std::wstring> SplitW(const std::wstring& s, wchar_t sep)
{
	std::vector<std::wstring> out;
	std::size_t start = 0;
	for (;;) {
		const std::size_t pos = s.find(sep, start);
		if (pos == std::wstring::npos) { out.push_back(s.substr(start)); break; }
		out.push_back(s.substr(start, pos - start));
		start = pos + 1;
	}
	return out;
}

std::wstring ToW(long long v)
{
	std::string s;
	AppendLong(s, v);
	return std::wstring(s.begin(), s.end());
}

//===========================================================================
// INI 섹션
//===========================================================================
struct RawSection {
	std::wstring title;
	std::size_t  headerLine;   // [제목] 줄 번호(0부터)
	std::size_t  firstLine;    // 데이터 첫 줄
	std::size_t  endLine;      // 다음 섹션 제목 줄(없으면 줄 수)
	std::vector<SXKeyValue> pairs;
};

void ParseSections(const std::vector<Line>& lines, std::vector<RawSection>& secs)
{
	for (std::size_t k = 0; k < lines.size(); ++k) {
		Line ln = lines[k];
		while (ln.e > ln.b && IsSpace(ln.e[-1])) --ln.e;
		const std::size_t n = static_cast<std::size_t>(ln.e - ln.b);
		if (n >= 2 && ln.b[0] == '[' && ln.e[-1] == ']') {
			if (!secs.empty()) secs.back().endLine = k;
			RawSection s;
			s.title = DecodeText(ln.b + 1, n - 2);
			s.headerLine = k;
			s.firstLine = k + 1;
			s.endLine = lines.size();
			secs.push_back(s);
			continue;
		}
		if (secs.empty() || n == 0 || ln.b[0] == ';') continue;
		const char* eq = static_cast<const char*>(std::memchr(ln.b, '=', n));
		if (!eq) continue;   // 데이터 줄
		SXKeyValue kv;
		kv.first = Trim(DecodeText(ln.b, static_cast<std::size_t>(eq - ln.b)));
		kv.second = DecodeText(eq + 1, static_cast<std::size_t>(ln.e - eq - 1));
		secs.back().pairs.push_back(kv);
	}
}

const std::wstring* FindValue(const RawSection& s, const wchar_t* key)
{
	for (std::size_t i = 0; i < s.pairs.size(); ++i)
		if (s.pairs[i].first == key) return &s.pairs[i].second;
	return 0;
}

bool WToFloat(const std::wstring& w, float& f)
{
	const std::string a = EncodeUtf8(Trim(w));
	return ParseFloat(a.data(), a.data() + a.size(), f);
}

bool WToLong(const std::wstring& w, long long& v)
{
	const std::string a = EncodeUtf8(Trim(w));
	return ParseLong(a.data(), a.data() + a.size(), v);
}

bool WToInt(const std::wstring& w, int& v)
{
	long long t;
	if (!WToLong(w, t) || t < INT_MIN || t > INT_MAX) return false;
	v = static_cast<int>(t);
	return true;
}

// "r,g,b,a" → 최대 4개. 모자라면 남는 칸은 그대로.
bool WToVec4(const std::wstring& w, float* v)
{
	std::vector<std::wstring> parts = SplitW(w, L',');
	if (!parts.empty() && Trim(parts.back()).empty()) parts.pop_back();   // 빈 값 · 끝 쉼표
	if (parts.size() > 4) return false;
	for (std::size_t i = 0; i < parts.size(); ++i)
		if (!WToFloat(parts[i], v[i])) return false;
	return true;
}

struct VarSpec {
	std::wstring name;
	int nRow;
	int nCol;
};

// "이름/행/열,이름/행/열"
bool ParseVarList(const std::wstring& value, std::vector<VarSpec>& out, std::wstring& err, const wchar_t* key)
{
	out.clear();
	if (Trim(value).empty()) return true;
	const std::vector<std::wstring> items = SplitW(value, L',');
	for (std::size_t i = 0; i < items.size(); ++i) {
		if (Trim(items[i]).empty() && i + 1 == items.size()) break;   // 끝 쉼표
		const std::vector<std::wstring> f = SplitW(items[i], L'/');
		VarSpec v;
		if (f.size() != 3 || !WToInt(f[1], v.nRow) || !WToInt(f[2], v.nCol) || v.nRow < 0 || v.nCol < 0) {
			err = std::wstring(L"[Info] ") + key + L" 항목 형식이 잘못됨: " + items[i];
			return false;
		}
		v.name = Trim(f[0]);
		out.push_back(v);
	}
	return true;
}

template <class T>
bool ReadRows(const std::vector<Line>& lines, const RawSection& sec, const VarSpec& spec, std::size_t bufSize,
              SXNamedContainer<T>& out, std::wstring& err)
{
	const std::size_t avail = sec.endLine - sec.firstLine;
	std::vector<Token> toks;

	// 옛 판의 1열 변수: 한 줄에 행 값을 모두 쉼표로 이어 썼다.
	if (spec.nCol == 1 && spec.nRow > 1 && avail >= 1) {
		const Line& ln = lines[sec.firstLine];
		if (!Tokenize(ln, toks)) { err = L"[" + spec.name + L"] 첫 줄에 빈 값이 있음"; return false; }
		if (toks.size() > 1) {
			if (toks.size() != static_cast<std::size_t>(spec.nRow)) {
				err = L"[" + spec.name + L"] 옛 1열 형식: 값 " + ToW(spec.nRow) + L"개가 필요한데 " +
				      ToW(static_cast<long long>(toks.size())) + L"개";
				return false;
			}
			out.Allocate(spec.nRow, 1);
			for (int r = 0; r < spec.nRow; ++r) {
				T v;
				if (!ParseValue(toks[r].b, toks[r].e, v)) {
					err = L"[" + spec.name + L"] " + ToW(r) + L"번째 값을 숫자로 읽을 수 없음";
					return false;
				}
				out[r][0] = v;
			}
			return true;
		}
	}

	if (avail < static_cast<std::size_t>(spec.nRow)) {
		err = L"[" + spec.name + L"] 행 " + ToW(spec.nRow) + L"개가 필요한데 " + ToW(static_cast<long long>(avail)) + L"줄뿐";
		return false;
	}

	// 1차: 행마다 값 개수 → 열 수 = max(헤더 열 수, 가장 긴 행)
	int maxCount = 0;
	for (int r = 0; r < spec.nRow; ++r) {
		if (!Tokenize(lines[sec.firstLine + r], toks)) {
			err = L"[" + spec.name + L"] " + ToW(static_cast<long long>(sec.firstLine + r + 1)) + L"번째 줄에 빈 값이 있음";
			return false;
		}
		if (static_cast<int>(toks.size()) > maxCount) maxCount = static_cast<int>(toks.size());
	}
	const int nCol = (spec.nCol > maxCount) ? spec.nCol : maxCount;
	// 헤더 열 수가 실제 행보다 터무니없이 크면(깨진 파일) 큰 메모리를 잡기 전에 멈춘다.
	if (nCol > maxCount && static_cast<double>(spec.nRow) * nCol > 16.0 * static_cast<double>(bufSize) + 4096.0) {
		err = L"[" + spec.name + L"] 헤더의 열 수(" + ToW(spec.nCol) + L")가 파일 내용에 비해 비정상";
		return false;
	}
	out.Allocate(spec.nRow, nCol);

	// 2차: 값 읽기
	for (int r = 0; r < spec.nRow; ++r) {
		Tokenize(lines[sec.firstLine + r], toks);
		SXRow<T> row = out[r];
		for (std::size_t j = 0; j < toks.size(); ++j) {
			if (!ParseValue(toks[j].b, toks[j].e, row.rpData[j])) {
				err = L"[" + spec.name + L"] " + ToW(static_cast<long long>(sec.firstLine + r + 1)) +
				      L"번째 줄의 값을 숫자로 읽을 수 없음: " + DecodeText(toks[j].b, static_cast<std::size_t>(toks[j].e - toks[j].b));
				return false;
			}
		}
		out.SetRowSize(r, static_cast<int>(toks.size()));
	}
	return true;
}

bool StartsWith(const std::wstring& s, const wchar_t* prefix)
{
	const std::size_t n = std::wcslen(prefix);
	return s.size() >= n && s.compare(0, n, prefix) == 0;
}

bool HasLineBreak(const std::wstring& s)
{
	return s.find(L'\n') != std::wstring::npos || s.find(L'\r') != std::wstring::npos;
}

void AppendLine(std::string& out, const wchar_t* key, const std::wstring& value)
{
	out += EncodeUtf8(key);
	out += '=';
	out += EncodeUtf8(value);
	out += '\n';
}

template <class T>
void AppendRows(std::string& out, const SXNamedContainer<T>& c)
{
	out += '[';
	out += EncodeUtf8(c.Name);
	out += "]\n";
	for (int r = 0; r < c.size(); ++r) {
		SXRow<const T> row = c[r];
		for (int j = 0; j < row.npData; ++j) {
			if (j) out += ',';
			AppendValue(out, row.rpData[j]);
		}
		out += '\n';
	}
}

void AppendVec4(std::string& out, const char* key, const float* v)
{
	out += key;
	out += '=';
	for (int i = 0; i < 4; ++i) {
		if (i) out += ',';
		AppendFloat(out, v[i]);
	}
	out += '\n';
}

} // namespace

//===========================================================================
// SXMeshMaterial
//===========================================================================
SXMeshMaterial::SXMeshMaterial() : ID(0), Shininess(50.f), GLAlpha(1.f), tex_Kd(0), tex_Ks(0), tex_d(0), tex_bump(0)
{
	for (int i = 0; i < 4; ++i) { Ka[i] = 1.f; Kd[i] = 1.f; Ks[i] = 1.f; }   // 옛 SMTLData::Init 과 같은 기본값
}

//===========================================================================
// SXMesh — 변수
//===========================================================================
SXMesh::SXMesh() {}

bool SXMesh::IsValidVariableName(const std::wstring& name)
{
	if (name.empty()) return false;
	if (name[0] == L' ' || name[0] == L'\t' || name[name.size() - 1] == L' ' || name[name.size() - 1] == L'\t') return false;
	for (std::size_t i = 0; i < name.size(); ++i) {
		const wchar_t c = name[i];
		if (c < 0x20 || c == 0x7F || c == L'[' || c == L']' || c == L'=' || c == L',' || c == L'/') return false;
	}
	if (name == L"Info") return false;
	if (StartsWith(name, L"pMTLData(")) return false;
	return true;
}

void SXMesh::CheckNewName(const std::wstring& name, bool forNode) const
{
	if (!IsValidVariableName(name)) throw std::invalid_argument("SXMesh: invalid variable name");
	if (forNode ? (FindElementV(name) != 0) : (FindNodeV(name) != 0))
		throw std::invalid_argument("SXMesh: the name is already used by the other kind of variable");
	for (std::size_t i = 0; i < ExtraSections.size(); ++i)
		if (ExtraSections[i].Title == name) throw std::invalid_argument("SXMesh: the name is used by a preserved section");
}

SXNodeNamedContainer& SXMesh::NewNodeV(const std::wstring& name, int nRow, int nCol)
{
	if (SXNodeNamedContainer* c = FindNodeV(name)) { c->Allocate(nRow, nCol); return *c; }
	CheckNewName(name, true);
	SXNodeNamedContainer c(name);
	c.Allocate(nRow, nCol);
	NodeVs_.push_back(std::move(c));
	return NodeVs_.back();
}

SXElementNamedContainer& SXMesh::NewElementV(const std::wstring& name, int nRow, int nCol)
{
	if (SXElementNamedContainer* c = FindElementV(name)) { c->Allocate(nRow, nCol); return *c; }
	CheckNewName(name, false);
	SXElementNamedContainer c(name);
	c.Allocate(nRow, nCol);
	ElementVs_.push_back(std::move(c));
	return ElementVs_.back();
}

SXNodeNamedContainer* SXMesh::FindNodeV(const std::wstring& name)
{
	for (std::deque<SXNodeNamedContainer>::iterator it = NodeVs_.begin(); it != NodeVs_.end(); ++it)
		if (it->Name == name) return &*it;
	return 0;
}

const SXNodeNamedContainer* SXMesh::FindNodeV(const std::wstring& name) const
{
	for (std::deque<SXNodeNamedContainer>::const_iterator it = NodeVs_.begin(); it != NodeVs_.end(); ++it)
		if (it->Name == name) return &*it;
	return 0;
}

SXElementNamedContainer* SXMesh::FindElementV(const std::wstring& name)
{
	for (std::deque<SXElementNamedContainer>::iterator it = ElementVs_.begin(); it != ElementVs_.end(); ++it)
		if (it->Name == name) return &*it;
	return 0;
}

const SXElementNamedContainer* SXMesh::FindElementV(const std::wstring& name) const
{
	for (std::deque<SXElementNamedContainer>::const_iterator it = ElementVs_.begin(); it != ElementVs_.end(); ++it)
		if (it->Name == name) return &*it;
	return 0;
}

SXNodeNamedContainer& SXMesh::GetNodeV(const std::wstring& name)
{
	if (SXNodeNamedContainer* c = FindNodeV(name)) return *c;
	throw std::out_of_range("SXMesh::GetNodeV: no node variable with that name");
}

const SXNodeNamedContainer& SXMesh::GetNodeV(const std::wstring& name) const
{
	if (const SXNodeNamedContainer* c = FindNodeV(name)) return *c;
	throw std::out_of_range("SXMesh::GetNodeV: no node variable with that name");
}

SXElementNamedContainer& SXMesh::GetElementV(const std::wstring& name)
{
	if (SXElementNamedContainer* c = FindElementV(name)) return *c;
	throw std::out_of_range("SXMesh::GetElementV: no element variable with that name");
}

const SXElementNamedContainer& SXMesh::GetElementV(const std::wstring& name) const
{
	if (const SXElementNamedContainer* c = FindElementV(name)) return *c;
	throw std::out_of_range("SXMesh::GetElementV: no element variable with that name");
}

int SXMesh::GetNNode(const std::wstring& name) const
{
	const SXNodeNamedContainer* c = FindNodeV(name);
	return c ? c->size() : -1;
}

int SXMesh::GetNElement(const std::wstring& name) const
{
	const SXElementNamedContainer* c = FindElementV(name);
	return c ? c->size() : -1;
}

bool SXMesh::RemoveNodeV(const std::wstring& name)
{
	for (std::deque<SXNodeNamedContainer>::iterator it = NodeVs_.begin(); it != NodeVs_.end(); ++it)
		if (it->Name == name) { NodeVs_.erase(it); return true; }
	return false;
}

bool SXMesh::RemoveElementV(const std::wstring& name)
{
	for (std::deque<SXElementNamedContainer>::iterator it = ElementVs_.begin(); it != ElementVs_.end(); ++it)
		if (it->Name == name) { ElementVs_.erase(it); return true; }
	return false;
}

void SXMesh::Reset()
{
	NodeVs_.clear();
	ElementVs_.clear();
	Materials.clear();
	ExtraSections.clear();
	Display = SXMeshDisplay();
	Name.clear();
}

//===========================================================================
// SXMesh — 저장
//===========================================================================
bool SXMesh::SaveToStream(std::ostream& os, std::wstring* err) const
{
	std::wstring e;

	// 섹션 이름이 겹치거나 줄을 깨뜨리는 문자가 없는지 먼저 본다.
	std::set<std::wstring> titles;
	titles.insert(L"Info");
	for (std::size_t i = 0; i < NodeVs_.size() && e.empty(); ++i) {
		const std::wstring& n = NodeVs_[i].Name;
		if (!IsValidVariableName(n)) e = L"저장할 수 없는 변수 이름: " + n;
		else if (!titles.insert(n).second) e = L"변수 이름이 겹침: " + n;
	}
	for (std::size_t i = 0; i < ElementVs_.size() && e.empty(); ++i) {
		const std::wstring& n = ElementVs_[i].Name;
		if (!IsValidVariableName(n)) e = L"저장할 수 없는 변수 이름: " + n;
		else if (!titles.insert(n).second) e = L"변수 이름이 겹침: " + n;
	}
	for (std::size_t i = 0; i < ExtraSections.size() && e.empty(); ++i) {
		const SXIniSection& s = ExtraSections[i];
		if (s.Title.empty() || HasLineBreak(s.Title) || s.Title.find(L']') != std::wstring::npos || StartsWith(s.Title, L"pMTLData("))
			e = L"저장할 수 없는 섹션 이름: " + s.Title;
		else if (!titles.insert(s.Title).second)
			e = L"섹션 이름이 겹침: " + s.Title;
		for (std::size_t k = 0; k < s.Pairs.size() && e.empty(); ++k)
			if (s.Pairs[k].first.empty() || s.Pairs[k].first.find(L'=') != std::wstring::npos ||
			    HasLineBreak(s.Pairs[k].first) || HasLineBreak(s.Pairs[k].second))
				e = L"[" + s.Title + L"] 저장할 수 없는 키 또는 값: " + s.Pairs[k].first;
	}
	for (std::size_t m = 0; m < Materials.size() && e.empty(); ++m) {
		const SXMeshMaterial& t = Materials[m];
		if (HasLineBreak(t.Name) || HasLineBreak(t.map_Kd) || HasLineBreak(t.map_Ks) || HasLineBreak(t.map_d) || HasLineBreak(t.bump))
			e = L"재질 " + ToW(static_cast<long long>(m)) + L": 값에 줄바꿈이 있음";
		for (std::size_t k = 0; k < t.Extra.size() && e.empty(); ++k)
			if (t.Extra[k].first.empty() || t.Extra[k].first.find(L'=') != std::wstring::npos ||
			    HasLineBreak(t.Extra[k].first) || HasLineBreak(t.Extra[k].second))
				e = L"재질 " + ToW(static_cast<long long>(m)) + L": 저장할 수 없는 키 또는 값: " + t.Extra[k].first;
	}
	if (!e.empty()) { if (err) *err = e; return false; }

	std::string out;
	try {
		std::size_t estimate = 1024;
		for (std::size_t i = 0; i < NodeVs_.size(); ++i)    estimate += static_cast<std::size_t>(NodeVs_[i].size()) * static_cast<std::size_t>(NodeVs_[i].nCol()) * 12;
		for (std::size_t i = 0; i < ElementVs_.size(); ++i) estimate += static_cast<std::size_t>(ElementVs_[i].size()) * static_cast<std::size_t>(ElementVs_[i].nCol()) * 7;
		out.reserve(estimate);

		const SXNodeNamedContainer*    nodes = FindNodeV(L"ppNodes");
		const SXElementNamedContainer* elems = FindElementV(L"ppElements");
		const SXNodeNamedContainer*    uvw   = FindNodeV(L"ppTextureUVW");

		out += "[Info]\n";
		out += "FormatVersion=2\n";
		out += "nNode=";         AppendLong(out, nodes ? nodes->size() : 0);  out += '\n';
		out += "nElement=";      AppendLong(out, elems ? elems->size() : 0);  out += '\n';
		out += "nTexturePoint="; AppendLong(out, uvw ? uvw->size() : 0);      out += '\n';
		out += "iNodeShape=";    AppendLong(out, elems ? elems->nCol() : 0);  out += '\n';
		out += "GLAlpha=";       AppendFloat(out, Display.GLAlpha);           out += '\n';
		out += "nMTLData=";      AppendLong(out, static_cast<long long>(Materials.size())); out += '\n';
		out += "PropertyBit=";   AppendLong(out, static_cast<long long>(Display.PropertyBit)); out += '\n';
		out += "PointSize=";     AppendFloat(out, Display.PointSize);         out += '\n';
		out += "LineWidth=";     AppendFloat(out, Display.LineWidth);         out += '\n';

		out += "SXNNode=";
		for (std::size_t i = 0; i < NodeVs_.size(); ++i) {
			if (i) out += ',';
			out += EncodeUtf8(NodeVs_[i].Name);
			out += '/'; AppendLong(out, NodeVs_[i].size());
			out += '/'; AppendLong(out, NodeVs_[i].nCol());
		}
		out += '\n';
		out += "SXNElement=";
		for (std::size_t i = 0; i < ElementVs_.size(); ++i) {
			if (i) out += ',';
			out += EncodeUtf8(ElementVs_[i].Name);
			out += '/'; AppendLong(out, ElementVs_[i].size());
			out += '/'; AppendLong(out, ElementVs_[i].nCol());
		}
		out += '\n';

		for (std::size_t m = 0; m < Materials.size(); ++m) {
			const SXMeshMaterial& t = Materials[m];
			out += "[pMTLData("; AppendLong(out, static_cast<long long>(m)); out += ")]\n";
			out += "ID="; AppendLong(out, t.ID); out += '\n';
			AppendLine(out, L"Name", t.Name);
			AppendLine(out, L"map_Kd", t.map_Kd);
			AppendLine(out, L"map_Ks", t.map_Ks);
			AppendLine(out, L"map_d", t.map_d);
			AppendLine(out, L"bump", t.bump);
			AppendVec4(out, "Ka", t.Ka);
			AppendVec4(out, "Kd", t.Kd);
			AppendVec4(out, "Ks", t.Ks);
			out += "Shininess="; AppendFloat(out, t.Shininess); out += '\n';
			out += "GLAlpha=";   AppendFloat(out, t.GLAlpha);   out += '\n';
			for (std::size_t k = 0; k < t.Extra.size(); ++k) AppendLine(out, t.Extra[k].first.c_str(), t.Extra[k].second);
		}

		for (std::size_t i = 0; i < ExtraSections.size(); ++i) {
			const SXIniSection& s = ExtraSections[i];
			out += '['; out += EncodeUtf8(s.Title); out += "]\n";
			for (std::size_t k = 0; k < s.Pairs.size(); ++k) AppendLine(out, s.Pairs[k].first.c_str(), s.Pairs[k].second);
		}

		for (std::size_t i = 0; i < NodeVs_.size(); ++i)    AppendRows(out, NodeVs_[i]);
		for (std::size_t i = 0; i < ElementVs_.size(); ++i) AppendRows(out, ElementVs_[i]);
	} catch (const std::bad_alloc&) {
		if (err) *err = L"메모리가 모자라 저장하지 못함";
		return false;
	}

	os.write(out.data(), static_cast<std::streamsize>(out.size()));
	if (!os) { if (err) *err = L"스트림에 쓰지 못함"; return false; }
	return true;
}

bool SXMesh::SaveToFile(const std::wstring& path, std::wstring* err) const
{
#ifdef _MSC_VER
	std::ofstream f(path.c_str(), std::ios::binary | std::ios::trunc);
#else
	std::ofstream f(EncodeUtf8(path).c_str(), std::ios::binary | std::ios::trunc);
#endif
	if (!f) { if (err) *err = L"파일을 열 수 없음: " + path; return false; }
	if (!SaveToStream(f, err)) return false;
	f.close();
	if (!f) { if (err) *err = L"파일을 닫지 못함: " + path; return false; }
	return true;
}

//===========================================================================
// SXMesh — 불러오기
//===========================================================================
bool SXMesh::LoadFromBuffer(const char* buffer, std::size_t size, std::wstring* err)
{
	std::wstring e;
	try {
		if (buffer == 0 && size > 0) { if (err) *err = L"버퍼가 NULL"; return false; }

		std::vector<Line> lines;
		SplitLines(buffer, size, lines);
		std::vector<RawSection> secs;
		ParseSections(lines, secs);

		std::map<std::wstring, std::size_t> byTitle;   // 같은 제목이 여럿이면 처음 것
		for (std::size_t i = 0; i < secs.size(); ++i) byTitle.insert(std::make_pair(secs[i].title, i));
		std::vector<bool> used(secs.size(), false);

		std::map<std::wstring, std::size_t>::const_iterator infoIt = byTitle.find(L"Info");
		if (infoIt == byTitle.end()) { if (err) *err = L"[Info] 섹션이 없음(NMF 파일이 아님)"; return false; }
		const RawSection& info = secs[infoIt->second];
		used[infoIt->second] = true;

		SXMesh tmp;

		// 표시 옵션
		if (const std::wstring* v = FindValue(info, L"PropertyBit")) {
			long long t;
			if (!WToLong(*v, t) || t < 0 || t > 0xFFFFFFFFLL) { if (err) *err = L"[Info] PropertyBit 값이 잘못됨: " + *v; return false; }
			tmp.Display.PropertyBit = static_cast<unsigned long>(t);
		} else if (const std::wstring* b = FindValue(info, L"bData")) {
			long long t;   // 옛 판: bData=<PropertyBit 정수>
			if (WToLong(*b, t) && t >= 0 && t <= 0xFFFFFFFFLL) tmp.Display.PropertyBit = static_cast<unsigned long>(t);
		}
		if (const std::wstring* v = FindValue(info, L"PointSize")) WToFloat(*v, tmp.Display.PointSize);
		if (const std::wstring* v = FindValue(info, L"LineWidth")) WToFloat(*v, tmp.Display.LineWidth);
		if (const std::wstring* v = FindValue(info, L"GLAlpha"))   WToFloat(*v, tmp.Display.GLAlpha);

		// 변수
		std::vector<VarSpec> nodeSpecs, elemSpecs;
		if (const std::wstring* v = FindValue(info, L"SXNNode"))
			if (!ParseVarList(*v, nodeSpecs, e, L"SXNNode")) { if (err) *err = e; return false; }
		if (const std::wstring* v = FindValue(info, L"SXNElement"))
			if (!ParseVarList(*v, elemSpecs, e, L"SXNElement")) { if (err) *err = e; return false; }

		for (int pass = 0; pass < 2; ++pass) {
			const std::vector<VarSpec>& specs = pass == 0 ? nodeSpecs : elemSpecs;
			for (std::size_t i = 0; i < specs.size(); ++i) {
				const VarSpec& s = specs[i];
				if (!IsValidVariableName(s.name)) { if (err) *err = L"쓸 수 없는 변수 이름: " + s.name; return false; }
				if (tmp.FindNodeV(s.name) || tmp.FindElementV(s.name)) { if (err) *err = L"변수 이름이 겹침: " + s.name; return false; }
				std::map<std::wstring, std::size_t>::const_iterator it = byTitle.find(s.name);
				if (it == byTitle.end()) { if (err) *err = L"[" + s.name + L"] 섹션이 없음"; return false; }
				used[it->second] = true;
				if (pass == 0) {
					SXNodeNamedContainer& c = tmp.NewNodeV(s.name, 0, 0);
					if (!ReadRows(lines, secs[it->second], s, size, c, e)) { if (err) *err = e; return false; }
				} else {
					SXElementNamedContainer& c = tmp.NewElementV(s.name, 0, 0);
					if (!ReadRows(lines, secs[it->second], s, size, c, e)) { if (err) *err = e; return false; }
				}
			}
		}

		// 재질
		std::vector<std::size_t> mtlSecs;   // pMTLData(...) 섹션, 파일 순서
		for (std::size_t i = 0; i < secs.size(); ++i)
			if (StartsWith(secs[i].title, L"pMTLData(") && byTitle[secs[i].title] == i) mtlSecs.push_back(i);
		int nMTL = static_cast<int>(mtlSecs.size());
		if (const std::wstring* v = FindValue(info, L"nMTLData"))
			if (!WToInt(*v, nMTL) || nMTL < 0) { if (err) *err = L"[Info] nMTLData 값이 잘못됨: " + *v; return false; }
		for (int m = 0; m < nMTL; ++m) {
			std::size_t si = secs.size();
			std::map<std::wstring, std::size_t>::const_iterator it = byTitle.find(L"pMTLData(" + ToW(m) + L")");
			if (it != byTitle.end() && !used[it->second]) si = it->second;
			for (std::size_t k = 0; k < mtlSecs.size() && si == secs.size(); ++k)   // 옛 판: 제목이 pMTLData(ID)
				if (!used[mtlSecs[k]]) si = mtlSecs[k];
			if (si == secs.size()) { if (err) *err = L"재질 " + ToW(m) + L" 섹션이 없음(nMTLData=" + ToW(nMTL) + L")"; return false; }
			used[si] = true;

			const RawSection& s = secs[si];
			SXMeshMaterial t;
			t.ID = m;
			{
				const std::wstring inner = s.title.substr(9, s.title.size() >= 10 ? s.title.size() - 10 : 0);   // "pMTLData(" … ")"
				int id;
				if (WToInt(inner, id)) t.ID = id;
			}
			for (std::size_t k = 0; k < s.pairs.size(); ++k) {
				const std::wstring& key = s.pairs[k].first;
				const std::wstring& val = s.pairs[k].second;
				bool ok = true;
				if (key == L"ID")             ok = WToInt(val, t.ID);
				else if (key == L"Name")      t.Name = val;
				else if (key == L"map_Kd")    t.map_Kd = val;
				else if (key == L"map_Ks")    t.map_Ks = val;
				else if (key == L"map_d")     t.map_d = val;
				else if (key == L"bump")      t.bump = val;
				else if (key == L"Ka")        ok = WToVec4(val, t.Ka);
				else if (key == L"Kd")        ok = WToVec4(val, t.Kd);
				else if (key == L"Ks")        ok = WToVec4(val, t.Ks);
				else if (key == L"Shininess") ok = WToFloat(val, t.Shininess);
				else if (key == L"GLAlpha")   ok = WToFloat(val, t.GLAlpha);
				else                          t.Extra.push_back(s.pairs[k]);
				if (!ok) { if (err) *err = L"[" + s.title + L"] " + key + L" 값이 잘못됨: " + val; return false; }
			}
			tmp.Materials.push_back(t);
		}

		// 나머지 섹션은 그대로 보존(예: [ManikinOBB])
		for (std::size_t i = 0; i < secs.size(); ++i) {
			if (used[i] || byTitle[secs[i].title] != i) continue;
			SXIniSection x;
			x.Title = secs[i].title;
			x.Pairs = secs[i].pairs;
			tmp.ExtraSections.push_back(x);
		}

		tmp.Name = Name;
		*this = std::move(tmp);
		return true;
	} catch (const std::bad_alloc&) {
		e = L"메모리가 모자라 불러오지 못함(행 · 열 수가 비정상일 수 있음)";
	} catch (const std::exception& ex) {
		const char* w = ex.what();
		e = L"불러오기 실패: " + DecodeText(w, std::strlen(w));
	}
	if (err) *err = e;
	return false;
}

bool SXMesh::LoadFromStream(std::istream& is, std::wstring* err)
{
	std::string data;
	try {
		data.assign(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>());
	} catch (const std::exception&) {
		if (err) *err = L"스트림을 읽지 못함";
		return false;
	}
	if (is.bad()) { if (err) *err = L"스트림을 읽지 못함"; return false; }
	return LoadFromBuffer(data.data(), data.size(), err);
}

bool SXMesh::LoadFromFile(const std::wstring& path, std::wstring* err)
{
#ifdef _MSC_VER
	std::ifstream f(path.c_str(), std::ios::binary);
#else
	std::ifstream f(EncodeUtf8(path).c_str(), std::ios::binary);
#endif
	if (!f) { if (err) *err = L"파일을 열 수 없음: " + path; return false; }
	return LoadFromStream(f, err);
}
