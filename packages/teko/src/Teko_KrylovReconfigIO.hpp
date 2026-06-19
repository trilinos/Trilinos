// Teko_KrylovReconfigIO.hpp — the JSON file formats exchanged with the
// reconfiguration watcher (wait_for_request.py) and the per-solve result logs.
//
// Four file formats live here, all in the requests dir as s<N>_*.json:
//   * s<N>_request.json  (C++ -> watcher)  the surrogate C_hat / b_hat
//   * s<N>_reconfig.json (watcher -> C++)  the chosen ordering(s)
//   * s<N>_conv.json     (C++ output)      per-solve timings (solve1/solve2)
//   * s<N>_solved.json   (C++ output)      the test-ordering sweep
// plus the small record structs they serialize (ReconfigResponse, SolveStats,
// OrderingResult). There is no external JSON library — every file is written to
// a *.tmp and atomically renamed, and the reconfig response is read with a
// small targeted parser (not a general JSON reader). See ACTIVATION.md for the
// field-by-field documentation of each format.
//
// The reduced-model math (CHatData / computeCHat) lives in
// Teko_KrylovReducedModel.hpp, which this header includes for SDM / CHatData.
#pragma once

// ── Standard library ──────────────────────────────────────────────────────
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Filesystem (C++17)
#include <filesystem>
namespace fs = std::filesystem;

// ── Teuchos ───────────────────────────────────────────────────────────────
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

// ── Teko ──────────────────────────────────────────────────────────────────
#include "Teko_KrylovReducedModel.hpp"  // type aliases (SDM) + CHatData

namespace Teko {
namespace KrylovSurrogate {

// ═══════════════════════════════════════════════════════════════════════════
// s<N>_request.json — C++ → watcher (the surrogate)
// ═══════════════════════════════════════════════════════════════════════════

// Write a flat JSON array of ints: [a, b, c].
inline void writeIntArrayJson(std::ostream& os, const std::vector<int>& v)
{
    os << "[";
    for (std::size_t i = 0; i < v.size(); ++i) { if (i) os << ", "; os << v[i]; }
    os << "]";
}

// Write a dense matrix as a JSON array-of-arrays. (std::setprecision is sticky,
// so it is set once per stream — here and in the writers below.)
inline void writeDenseJson(std::ostream& os, const SDM& M, int indent)
{
    const std::string pad(indent, ' ');
    const std::string inner(indent + 2, ' ');
    os << std::setprecision(17) << pad << "[\n";
    for (int r = 0; r < M.numRows(); ++r) {
        os << inner << "[";
        for (int c = 0; c < M.numCols(); ++c) {
            if (c) os << ", ";
            os << M(r, c);
        }
        os << "]";
        if (r + 1 < M.numRows()) os << ",";
        os << "\n";
    }
    os << pad << "]";
}

// Write a (k × 1) column vector as a flat JSON array.
inline void writeVectorJson(std::ostream& os, const SDM& v)
{
    os << std::setprecision(17) << "[";
    for (int r = 0; r < v.numRows(); ++r) {
        if (r) os << ", ";
        os << v(r, 0);
    }
    os << "]";
}

// Scan requests_dir for existing s<N>_request.json files and return
// max(N) + 1, or 0 if none exist. This is how the C++ side picks the id for
// the next request, so request ids are unique and monotonically increasing
// across the lifetime of requests_dir (they are never reused, even if old
// s<N>_request.json files are deleted).
inline int nextRequestNumber(const std::string& requests_dir)
{
    int next = 0;
    if (!fs::exists(requests_dir)) return next;

    // Consider both s<N>_request.json and s<N>_conv.json. A non-converged
    // (record-only) solve writes a s<N>_conv.json with no matching
    // request.json; counting conv files too keeps that N reserved so a later
    // request can't reuse it and overwrite the record.
    const std::string prefix = "s";
    const std::array<std::string, 3> suffixes{"_request.json", "_conv.json", "_solved.json"};
    for (const auto& entry : fs::directory_iterator(requests_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string();
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        for (const auto& suffix : suffixes) {
            if (name.size() <= prefix.size() + suffix.size()) continue;
            if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0) continue;

            const std::string digits = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
            if (digits.empty() || !std::all_of(digits.begin(), digits.end(),
                                                [](char c){ return std::isdigit(static_cast<unsigned char>(c)); })) continue;

            const int n = std::stoi(digits);
            if (n + 1 > next) next = n + 1;
        }
    }
    return next;
}

// equation_ends[k] = sum(ranks[0..k-1]) for k = 0..nb, so equation_ends[0]==0
// and equation_ends[nb] == R (the total rank, R = sum(ranks)). Each
// equation_ends[k+1] is one-past-the-last row/col index of block k in the
// assembled C_hat matrix / b_hat vector below.
inline std::vector<int> computeEquationEnds(const CHatData& chat)
{
    std::vector<int> ends(chat.nb + 1, 0);
    for (int j = 0; j < chat.nb; ++j)
        ends[j + 1] = ends[j] + chat.ranks[j];
    return ends;
}

// Assemble the full R x R C_hat matrix (R = sum(ranks)) from chat.blocks,
// placing block (i,j) — size ranks[i] x ranks[j] — at rows
// [equation_ends[i], equation_ends[i+1]) and cols
// [equation_ends[j], equation_ends[j+1]).
inline SDM assembleCHat(const CHatData& chat, const std::vector<int>& equation_ends)
{
    const int R = equation_ends.back();
    SDM C(R, R);
    for (int i = 0; i < chat.nb; ++i) {
        for (int j = 0; j < chat.nb; ++j) {
            const SDM& blk = chat.blocks[i][j];
            for (int r = 0; r < blk.numRows(); ++r)
                for (int c = 0; c < blk.numCols(); ++c)
                    C(equation_ends[i] + r, equation_ends[j] + c) = blk(r, c);
        }
    }
    return C;
}

// Assemble the full R x 1 b_hat vector (R = sum(ranks)) by stacking
// chat.b_hat[j] — size ranks[j] x 1 — at rows
// [equation_ends[j], equation_ends[j+1]).
inline SDM assembleBHat(const CHatData& chat, const std::vector<int>& equation_ends)
{
    const int R = equation_ends.back();
    SDM b(R, 1);
    for (int j = 0; j < chat.nb; ++j) {
        const SDM& bj = chat.b_hat[j];
        for (int r = 0; r < bj.numRows(); ++r)
            b(equation_ends[j] + r, 0) = bj(r, 0);
    }
    return b;
}

// Write s<N>_request.json to requests_dir, where N = nextRequestNumber(requests_dir).
// Returns N so the caller can wait for the matching s<N>_reconfig.json
// and tag s<N>_conv.json with the same id.
//
// Format:
// {
//   "n_blocks": nb,
//   "block_sizes": [...],
//   "krylov_dim": m,
//   "ranks": [...],
//   "equation_ends": [0, r0, r0+r1, ..., R],   // length nb+1, R = sum(ranks)
//   "C_hat": [ [ ... ], ... ],                  // single R x R matrix
//   "b_hat": [ ... ]                            // single length-R vector
// }
//
// C_hat is the full R x R surrogate operator (R = sum(ranks)) and b_hat the
// full length-R reduced RHS, both assembled block-by-block. "equation_ends"
// gives the block boundaries within them: equation_ends[k] is the first
// row/col index of block k and equation_ends[k+1] one past its last, so
// block k occupies [equation_ends[k], equation_ends[k+1]).
//
// Every entry scales with nb (number of blocks) and ranks[j] (Krylov-derived
// rank, <= krylov_dim) only — nothing here scales with the global problem
// size n, so this stays small even for very large systems.
inline int writeRequestJson(const CHatData& chat, const std::string& requests_dir)
{
    const int request_id = nextRequestNumber(requests_dir);

    const std::vector<int> equation_ends = computeEquationEnds(chat);
    const SDM C = assembleCHat(chat, equation_ends);
    const SDM b = assembleBHat(chat, equation_ends);

    // Write to a temp file and atomically rename into place. Without this, a
    // reader's fs::exists(path) check (in wait_for_request.py) can observe
    // the file mid-write and parse incomplete JSON. A same-directory
    // rename() is atomic on POSIX, so readers only ever see the file fully
    // written or not at all.
    const std::string path     = requests_dir + "/s" + std::to_string(request_id) + "_request.json";
    const std::string tmp_path = path + ".tmp";
    std::ofstream ofs(tmp_path);
    TEUCHOS_TEST_FOR_EXCEPTION(!ofs.is_open(), std::runtime_error,
        "Teko::KrylovSurrogate::writeRequestJson: cannot open " + tmp_path);

    ofs << "{\n";
    ofs << "  \"n_blocks\": " << chat.nb << ",\n";
    ofs << "  \"block_sizes\": ";   writeIntArrayJson(ofs, chat.block_sizes); ofs << ",\n";
    ofs << "  \"krylov_dim\": " << chat.krylov_dim << ",\n";
    ofs << "  \"ranks\": ";         writeIntArrayJson(ofs, chat.ranks);       ofs << ",\n";
    ofs << "  \"equation_ends\": "; writeIntArrayJson(ofs, equation_ends);    ofs << ",\n";

    // C_hat: single R x R matrix
    ofs << "  \"C_hat\":\n";
    writeDenseJson(ofs, C, 2);
    ofs << ",\n";

    // b_hat: single length-R vector
    ofs << "  \"b_hat\": ";
    writeVectorJson(ofs, b);
    ofs << "\n";

    ofs << "}\n";
    ofs.close();

    fs::rename(tmp_path, path);

    std::cout << "[TekoAdaptive] wrote " << path << "\n";
    return request_id;
}

// ═══════════════════════════════════════════════════════════════════════════
// s<N>_reconfig.json — watcher → C++ (the response)
// ═══════════════════════════════════════════════════════════════════════════

// ── Minimal strict parsing for the reconfig response ──────────────────────
// s<N>_reconfig.json is machine-generated by wait_for_request.py in a fixed
// shape, e.g.
//   {"selection_mode": "chosen",
//    "use_ordering": [0,0,1],
//    "opt_ordering": [0,0,1],
//    "exh_opt_ordering": [],
//    "test_orderings": [[0,1,2],[0,0,1],[0,0,0]]}
// These helpers extract exactly the fields the C++ side consumes — a string
// enum, an int array, and an array-of-int-arrays — keyed off field names and
// rejecting malformed numeric tokens. It is deliberately a small targeted
// parser, not a general JSON reader, but keying off the field names makes it
// tolerant of extra fields in any position.

// Parse one integer from a token, rejecting anything std::stoi would silently
// truncate (e.g. "0.9" -> 0).
inline int parseStrictInt(std::string tok, const std::string& ctx)
{
    tok.erase(std::remove_if(tok.begin(), tok.end(),
                             [](char c){ return std::isspace(static_cast<unsigned char>(c)); }),
              tok.end());
    std::size_t consumed = 0;
    int value = 0;
    try {
        value = std::stoi(tok, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(ctx + ": non-integer entry \"" + tok + "\"");
    }
    if (consumed != tok.size())
        throw std::runtime_error(ctx + ": non-integer entry \"" + tok + "\"");
    return value;
}

// Index of the ']' matching the '[' at openPos (balanced). npos if unbalanced.
inline std::size_t matchBracket(const std::string& s, std::size_t openPos)
{
    int depth = 0;
    for (std::size_t i = openPos; i < s.size(); ++i) {
        if (s[i] == '[') ++depth;
        else if (s[i] == ']') { if (--depth == 0) return i; }
    }
    return std::string::npos;
}

// Parse "a, b, c" (comma-separated, no brackets) into ints.
inline std::vector<int> parseIntCsv(const std::string& body, const std::string& ctx)
{
    std::vector<int> out;
    std::istringstream ss(body);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        const bool nonspace = std::any_of(tok.begin(), tok.end(),
            [](char c){ return !std::isspace(static_cast<unsigned char>(c)); });
        if (nonspace) out.push_back(parseStrictInt(tok, ctx));
    }
    return out;
}

// Optional "key": "value" string field. Returns false if the key is absent.
inline bool parseStringField(const std::string& content, const std::string& key,
                             std::string& out)
{
    const auto kp = content.find("\"" + key + "\"");
    if (kp == std::string::npos) return false;
    const auto colon = content.find(':', kp + key.size() + 2);
    if (colon == std::string::npos) return false;
    const auto q1 = content.find('"', colon);
    if (q1 == std::string::npos) return false;
    const auto q2 = content.find('"', q1 + 1);
    if (q2 == std::string::npos) return false;
    out = content.substr(q1 + 1, q2 - q1 - 1);
    return true;
}

// Optional "key": [ints] field. Returns false if the key is absent.
inline bool parseIntArrayField(const std::string& content, const std::string& key,
                               std::vector<int>& out, const std::string& ctx)
{
    const auto kp = content.find("\"" + key + "\"");
    if (kp == std::string::npos) return false;
    const auto lb = content.find('[', kp);
    if (lb == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" has no array");
    const auto rb = matchBracket(content, lb);
    if (rb == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" array not closed");
    out = parseIntCsv(content.substr(lb + 1, rb - lb - 1), ctx + " \"" + key + "\"");
    return true;
}

// Optional "key": [[ints],[ints],...] field. Returns false if the key is absent.
inline bool parseListOfIntArraysField(const std::string& content, const std::string& key,
                                      std::vector<std::vector<int>>& out,
                                      const std::string& ctx)
{
    const auto kp = content.find("\"" + key + "\"");
    if (kp == std::string::npos) return false;
    const auto outerL = content.find('[', kp);
    if (outerL == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" has no array");
    const auto outerR = matchBracket(content, outerL);
    if (outerR == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" array not closed");
    out.clear();
    std::size_t i = outerL + 1;
    while (true) {
        const auto innerL = content.find('[', i);
        if (innerL == std::string::npos || innerL > outerR) break;
        const auto innerR = matchBracket(content, innerL);
        if (innerR == std::string::npos || innerR > outerR)
            throw std::runtime_error(ctx + ": \"" + key + "\" inner array not closed");
        out.push_back(parseIntCsv(content.substr(innerL + 1, innerR - innerL - 1),
                                  ctx + " \"" + key + "\""));
        i = innerR + 1;
    }
    return true;
}

// Parsed reconfig response. use_ordering empty => no response (timeout).
struct ReconfigResponse {
    std::string                   selection_mode;  // "chosen" | "best_conv" | "best_time"
    std::vector<int>              use_ordering;     // ordering applied in "chosen" mode
    std::vector<std::vector<int>> test_orderings;   // candidates to sweep/time (may be empty)
};

// Poll for s<request_id>_reconfig.json and parse it. The file is left in place
// after being read (NOT deleted) — its presence is how wait_for_request.py
// recognizes an already-answered request. Returns a response with empty
// use_ordering on timeout.
inline ReconfigResponse waitForReconfig(
    const std::string& requests_dir,
    int                request_id,
    int                timeout_s = 600)
{
    const std::string path = requests_dir + "/s" + std::to_string(request_id) + "_reconfig.json";
    const auto deadline =
        std::chrono::steady_clock::now() + std::chrono::seconds(timeout_s);

    std::cout << "[TekoAdaptive] waiting for " << path << " ...\n";
    std::cout.flush();

    while (std::chrono::steady_clock::now() < deadline) {
        if (fs::exists(path)) {
            std::ifstream ifs(path);
            std::string content((std::istreambuf_iterator<char>(ifs)),
                                 std::istreambuf_iterator<char>());
            ifs.close();

            const std::string ctx = "[TekoAdaptive] " + path;
            ReconfigResponse r;
            if (!parseIntArrayField(content, "use_ordering", r.use_ordering, ctx))
                throw std::runtime_error(ctx + ": no \"use_ordering\" field found");
            if (!parseStringField(content, "selection_mode", r.selection_mode))
                r.selection_mode = "chosen";
            parseListOfIntArraysField(content, "test_orderings", r.test_orderings, ctx);

            std::cout << "[TekoAdaptive] reconfig: mode=" << r.selection_mode
                      << "  use_ordering=[";
            for (int k = 0; k < (int)r.use_ordering.size(); ++k) {
                if (k) std::cout << ", ";
                std::cout << r.use_ordering[k];
            }
            std::cout << "]  test_orderings=" << r.test_orderings.size() << "\n";
            return r;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }

    std::cerr << "[TekoAdaptive] timed out waiting for " << path << "\n";
    return {};
}

// ═══════════════════════════════════════════════════════════════════════════
// s<N>_conv.json — C++ output (per-solve timings)
// ═══════════════════════════════════════════════════════════════════════════

// Convergence diagnostics for one FGMRES solve, written to
// s<N>_conv.json by writeConvergenceJson().
struct SolveStats {
    int    iterations;
    double initial_residual;
    double final_residual;
    double factor_wall_time_sec;
    double iterate_wall_time_sec;
    double total_wall_time_sec;
};

// Write s<request_id>_conv.json (atomically, via a .tmp file + rename).
// s2 == nullptr means no second solve was attempted (Phase 2 timed out waiting
// for s<request_id>_reconfig.json), and is recorded as "solve2": null.
inline void writeConvergenceJson(
    const std::string& convergence_dir,
    int                 request_id,
    const SolveStats&   s1,
    const SolveStats*   s2)
{
    fs::create_directories(convergence_dir);

    const std::string path     = convergence_dir + "/s" + std::to_string(request_id) + "_conv.json";
    const std::string tmp_path = path + ".tmp";
    std::ofstream ofs(tmp_path);
    TEUCHOS_TEST_FOR_EXCEPTION(!ofs.is_open(), std::runtime_error,
        "Teko::KrylovSurrogate::writeConvergenceJson: cannot open " + tmp_path);
    ofs << std::setprecision(17);

    auto writeStats = [&](const char* name, const SolveStats& s, bool trailing_comma) {
        ofs << "  \"" << name << "\": {\n";
        ofs << "    \"iterations\": " << s.iterations << ",\n";
        ofs << "    \"initial_residual\": " << s.initial_residual << ",\n";
        ofs << "    \"final_residual\": " << s.final_residual << ",\n";
        ofs << "    \"factor_wall_time_sec\": " << s.factor_wall_time_sec << ",\n";
        ofs << "    \"iterate_wall_time_sec\": " << s.iterate_wall_time_sec << ",\n";
        ofs << "    \"total_wall_time_sec\": " << s.total_wall_time_sec << ",\n";
        ofs << "    \"wall_time_sec\": " << s.total_wall_time_sec << "\n";
        ofs << "  }" << (trailing_comma ? ",\n" : "\n");
    };

    ofs << "{\n";
    ofs << "  \"request_id\": " << request_id << ",\n";
    writeStats("solve1", s1, /*trailing_comma=*/true);
    if (s2 != nullptr)
        writeStats("solve2", *s2, /*trailing_comma=*/false);
    else
        ofs << "  \"solve2\": null\n";
    ofs << "}\n";
    ofs.close();

    fs::rename(tmp_path, path);
    std::cout << "[TekoAdaptive] wrote " << path << "\n";
}

// ═══════════════════════════════════════════════════════════════════════════
// s<N>_solved.json — C++ output (the test-ordering sweep)
// ═══════════════════════════════════════════════════════════════════════════

// Per-ordering solve outcome (one row of s<N>_solved.json).
struct OrderingResult {
    std::vector<int> ordering;
    int    iterations           = 0;
    bool   converged            = false;
    double initial_residual     = 0.0;
    double final_residual       = 0.0;
    double factor_wall_time_sec = 0.0;
    double iterate_wall_time_sec = 0.0;
    double total_wall_time_sec  = 0.0;
};

// Write s<request_id>_solved.json (atomically): the full test-ordering sweep
// with per-ordering results, the selection_mode, and which ordering was
// selected/returned.
inline void writeSolvedJson(
    const std::string&                 dir,
    int                                request_id,
    const std::string&                 selection_mode,
    const std::vector<OrderingResult>& results,
    int                                selected_index,
    const std::vector<int>&            selected_ordering)
{
    fs::create_directories(dir);
    const std::string path = dir + "/s" + std::to_string(request_id) + "_solved.json";
    const std::string tmp_path = path + ".tmp";
    std::ofstream ofs(tmp_path);
    TEUCHOS_TEST_FOR_EXCEPTION(!ofs.is_open(), std::runtime_error,
        "Teko::KrylovSurrogate::writeSolvedJson: cannot open " + tmp_path);
    ofs << std::setprecision(17);

    ofs << "{\n";
    ofs << "  \"selection_mode\": \"" << selection_mode << "\",\n";
    ofs << "  \"request_id\": " << request_id << ",\n";
    ofs << "  \"selected_index\": " << selected_index << ",\n";
    ofs << "  \"selected_ordering\": "; writeIntArrayJson(ofs, selected_ordering); ofs << ",\n";
    ofs << "  \"results\": [\n";
    for (std::size_t i = 0; i < results.size(); ++i) {
        const OrderingResult& r = results[i];
        ofs << "    {\n";
        ofs << "      \"ordering\": "; writeIntArrayJson(ofs, r.ordering); ofs << ",\n";
        ofs << "      \"iterations\": " << r.iterations << ",\n";
        ofs << "      \"converged\": " << (r.converged ? "true" : "false") << ",\n";
        ofs << "      \"initial_residual\": " << r.initial_residual << ",\n";
        ofs << "      \"final_residual\": " << r.final_residual << ",\n";
        ofs << "      \"factor_wall_time_sec\": " << r.factor_wall_time_sec << ",\n";
        ofs << "      \"iterate_wall_time_sec\": " << r.iterate_wall_time_sec << ",\n";
        ofs << "      \"total_wall_time_sec\": " << r.total_wall_time_sec << "\n";
        ofs << "    }" << (i + 1 < results.size() ? "," : "") << "\n";
    }
    ofs << "  ]\n";
    ofs << "}\n";
    ofs.close();

    fs::rename(tmp_path, path);
    std::cout << "[TekoAdaptive] wrote " << path << "\n";
}

} // namespace KrylovSurrogate
} // namespace Teko
