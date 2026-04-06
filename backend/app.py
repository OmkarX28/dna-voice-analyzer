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

if __name__ == '__main__':
    app.run(debug=True)