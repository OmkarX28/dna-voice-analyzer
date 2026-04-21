from flask import Flask, request, jsonify
from flask_cors import CORS
from dotenv import load_dotenv
from Bio.Seq import Seq
import os

load_dotenv()

app = Flask(__name__)
CORS(app)

# ─── KMP Algorithm ───────────────────────────────────────────
def build_lps(pattern):
    lps = [0] * len(pattern)
    length = 0
    i = 1
    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def kmp_search(sequence, pattern):
    positions = []
    n, m = len(sequence), len(pattern)
    lps = build_lps(pattern)
    i = j = 0
    while i < n:
        if sequence[i] == pattern[j]:
            i += 1
            j += 1
        if j == m:
            positions.append(i - j)
            j = lps[j - 1]
        elif i < n and sequence[i] != pattern[j]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return positions

# ─── GC Content ──────────────────────────────────────────────
def calculate_gc(sequence):
    g = sequence.count('G')
    c = sequence.count('C')
    total = len(sequence)
    if total == 0:
        return 0
    return round((g + c) / total * 100, 2)

# ─── Routes ──────────────────────────────────────────────────
@app.route('/ping', methods=['GET'])
def ping():
    return jsonify({"message": "DNA Analyzer backend is running!"})

@app.route('/analyze/pattern', methods=['POST'])
def find_pattern():
    data = request.get_json()
    sequence = data.get('sequence', '').upper().strip()
    pattern = data.get('pattern', '').upper().strip()

    if not sequence or not pattern:
        return jsonify({"error": "Please provide both a sequence and a pattern"}), 400

    if not all(c in 'ATCG' for c in sequence):
        return jsonify({"error": "Sequence contains invalid characters. Only A, T, C, G allowed"}), 400

    if not all(c in 'ATCG' for c in pattern):
        return jsonify({"error": "Pattern contains invalid characters. Only A, T, C, G allowed"}), 400

    positions = kmp_search(sequence, pattern)
    gc = calculate_gc(sequence)

    return jsonify({
        "sequence_length": len(sequence),
        "pattern": pattern,
        "pattern_length": len(pattern),
        "occurrences": len(positions),
        "positions": positions,
        "gc_content": gc,
        "message": f"Pattern '{pattern}' found {len(positions)} time(s) in the sequence." if positions else f"Pattern '{pattern}' was not found in the sequence."
    })
# ─── Mutation Detector ───────────────────────────────────────
def detect_mutations(reference, sample):
    mutations = []

    # Point mutations (substitutions)
    min_len = min(len(reference), len(sample))
    for i in range(min_len):
        if reference[i] != sample[i]:
            mutations.append({
                "type": "point_mutation",
                "position": i,
                "reference_base": reference[i],
                "sample_base": sample[i],
                "description": f"Position {i}: '{reference[i]}' changed to '{sample[i]}'"
            })

    # Insertion (sample is longer than reference)
    if len(sample) > len(reference):
        inserted = sample[len(reference):]
        mutations.append({
            "type": "insertion",
            "position": len(reference),
            "inserted_bases": inserted,
            "description": f"Insertion of '{inserted}' after position {len(reference) - 1}"
        })

    # Deletion (sample is shorter than reference)
    if len(sample) < len(reference):
        deleted = reference[len(sample):]
        mutations.append({
            "type": "deletion",
            "position": len(sample),
            "deleted_bases": deleted,
            "description": f"Deletion of '{deleted}' starting at position {len(sample)}"
        })

    return mutations


@app.route('/analyze/mutations', methods=['POST'])
def find_mutations():
    data = request.get_json()
    reference = data.get('reference', '').upper().strip()
    sample = data.get('sample', '').upper().strip()

    if not reference or not sample:
        return jsonify({"error": "Please provide both a reference and a sample sequence"}), 400

    if not all(c in 'ATCG' for c in reference):
        return jsonify({"error": "Reference contains invalid characters. Only A, T, C, G allowed"}), 400

    if not all(c in 'ATCG' for c in sample):
        return jsonify({"error": "Sample contains invalid characters. Only A, T, C, G allowed"}), 400

    mutations = detect_mutations(reference, sample)
    gc_ref = calculate_gc(reference)
    gc_sample = calculate_gc(sample)

    return jsonify({
        "reference_length": len(reference),
        "sample_length": len(sample),
        "total_mutations": len(mutations),
        "mutations": mutations,
        "gc_content": {
            "reference": gc_ref,
            "sample": gc_sample
        },
        "message": f"{len(mutations)} mutation(s) detected." if mutations else "No mutations detected. Sequences are identical."
    })
# ─── Smith-Waterman Alignment ────────────────────────────────
def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    matrix = [[0] * cols for _ in range(rows)]
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diagonal = matrix[i-1][j-1] + match_score
            up = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap
            matrix[i][j] = max(0, diagonal, up, left)
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_pos

    while i > 0 and j > 0 and matrix[i][j] > 0:
        current = matrix[i][j]
        match_score = match if seq1[i-1] == seq2[j-1] else mismatch
        if current == matrix[i-1][j-1] + match_score:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif current == matrix[i-1][j] + gap:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    similarity = round((matches / max(len(aligned_seq1), 1)) * 100, 2)

    return max_score, aligned_seq1, aligned_seq2, similarity


@app.route('/analyze/align', methods=['POST'])
def align_sequences():
    data = request.get_json()
    seq1 = data.get('seq1', '').upper().strip()
    seq2 = data.get('seq2', '').upper().strip()

    if not seq1 or not seq2:
        return jsonify({"error": "Please provide both seq1 and seq2"}), 400

    if not all(c in 'ATCG' for c in seq1):
        return jsonify({"error": "seq1 contains invalid characters. Only A, T, C, G allowed"}), 400

    if not all(c in 'ATCG' for c in seq2):
        return jsonify({"error": "seq2 contains invalid characters. Only A, T, C, G allowed"}), 400

    score, aligned1, aligned2, similarity = smith_waterman(seq1, seq2)

    return jsonify({
        "seq1_length": len(seq1),
        "seq2_length": len(seq2),
        "alignment_score": score,
        "similarity_percentage": similarity,
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
        "message": f"Sequences are {similarity}% similar with an alignment score of {score}."
    })

if __name__ == '__main__':
    app.run(debug=True)