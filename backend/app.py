from flask import Flask, request, jsonify, make_response
from flask_cors import CORS
from dotenv import load_dotenv
from Bio.Seq import Seq
from groq import Groq
import json
import os

from pathlib import Path
load_dotenv(dotenv_path=Path(__file__).resolve().parent.parent / ".env", override=True)
# ── Groq Setup ─────────────────────────────────────────────
groq_client = Groq(api_key=os.getenv("GROQ_API_KEY"))

app = Flask(__name__)
CORS(app)
@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type')
    response.headers.add('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
    return response
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

    # Build scoring matrix
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

    # Traceback to find aligned sequences
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

    # Calculate similarity percentage
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
# ── NLP Intent Engine ─────────────────────────────────────────
def parse_intent(user_query, sequence=None):
    prompt = f"""
You are an AI assistant for a DNA sequence analysis tool.
A user has sent the following voice query: "{user_query}"
The DNA sequence they uploaded is: "{sequence if sequence else 'not provided'}"

Your job is to:
1. Understand what analysis the user wants
2. Extract any pattern or second sequence they mentioned
3. Return ONLY a valid JSON object — no explanation, no markdown, no extra text

Return this exact JSON structure:
{{
  "intent": "<one of: pattern_search, mutation_detection, alignment, gc_content, explain, unknown>",
  "pattern": "<the DNA pattern to search for, or null>",
  "second_sequence": "<a second sequence for alignment or mutation comparison, or null>",
  "plain_english_query": "<rephrase what the user wants in one clear sentence>"
}}
"""
    response = groq_client.chat.completions.create(
        model="llama-3.3-70b-versatile",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.1
    )
    text = response.choices[0].message.content.strip()
    text = text.replace("```json", "").replace("```", "").strip()
    return json.loads(text)


MUTATION_KNOWLEDGE_BASE = """
KNOWN MUTATIONS REFERENCE:
1. Sickle Cell Anemia - HBB gene - Point mutation GAG→GTG at codon 6
   Causes: Red blood cells become sickle-shaped, blocking blood flow
   Explanation: The single change from GAG to GTG causes hemoglobin to form rigid rods inside red blood cells, making them crescent-shaped instead of round. This causes severe pain, anemia, and organ damage.
   Severity: High

2. BRCA1 Breast Cancer - BRCA1 gene - Deletion/frameshift mutations
   Causes: Significantly increased risk of breast and ovarian cancer
   Explanation: Mutations in BRCA1 disable a tumor suppressor protein that normally repairs damaged DNA. Without this repair mechanism, cells can become cancerous. Carriers have up to 72% lifetime risk of breast cancer.
   Severity: High

3. Cystic Fibrosis - CFTR gene - Deletion of phenylalanine at position 508
   Causes: Thick mucus buildup in lungs, digestive system, and other organs
   Explanation: The missing amino acid causes the CFTR protein to fold incorrectly and get destroyed before reaching the cell surface. Without functional CFTR, chloride transport fails causing thick sticky mucus to accumulate.
   Severity: High

4. Huntington's Disease - HTT gene - CAG repeat expansion beyond 36 repeats
   Causes: Progressive brain cell death leading to movement, cognitive and psychiatric disorders
   Explanation: Extra CAG repeats create an abnormally long huntingtin protein that damages brain cells over time. Symptoms typically appear between ages 30-50 and worsen progressively with no cure currently available.
   Severity: High

5. TB Drug Resistance - rpoB gene - Point mutations at codons 516, 526, 531
   Causes: Resistance to rifampicin, a key TB antibiotic
   Explanation: Mutations in the rpoB gene change the shape of RNA polymerase so rifampicin can no longer bind to it. This makes the TB bacteria resistant to one of the most important first-line antibiotics, requiring more toxic second-line drugs.
   Severity: High

6. Beta Thalassemia - HBB gene - Various point mutations and deletions
   Causes: Reduced or absent beta-globin production leading to severe anemia
   Explanation: Mutations in the HBB gene reduce production of the beta-globin component of hemoglobin. Without enough functional hemoglobin, red blood cells are small and destroyed rapidly, causing severe anemia requiring regular blood transfusions.
   Severity: High

7. PKU (Phenylketonuria) - PAH gene - Point mutations reducing enzyme activity
   Causes: Buildup of phenylalanine in blood causing intellectual disability
   Explanation: Mutations in the PAH gene reduce or eliminate the enzyme that breaks down phenylalanine from food. Without treatment, phenylalanine accumulates to toxic levels in the brain causing severe intellectual disability. Treatable with a strict low-phenylalanine diet.
   Severity: High

8. Hemophilia A - F8 gene - Inversion or point mutations
   Causes: Inability to form blood clots properly
   Explanation: Mutations in the F8 gene reduce or eliminate clotting factor VIII. Without this protein, the blood clotting cascade cannot complete properly, causing prolonged bleeding from injuries and spontaneous bleeding into joints and muscles.
   Severity: High
"""

def generate_plain_english_response(intent, analysis_result, user_query, context=''):
    prompt = f"""
You are a helpful medical AI assistant explaining DNA analysis results.
You have access to this medical knowledge base about known mutations:

{MUTATION_KNOWLEDGE_BASE}

The user asked: "{user_query}"
The analysis result is: {json.dumps(analysis_result)}
IMPORTANT CONTEXT FROM THE ANALYSIS: {context}

Based on the context above, the analysis HAS already been run and found specific results.
Do NOT say no mutations were found if the context says mutations were found.
Do NOT say you need more information — you already have it in the context.

Write a clear, friendly, plain-English explanation in 2-3 sentences.
If the context mentions 1 mutation was found and the sequence matches the sickle cell mutation pattern, explain it specifically.
Speak directly to the user. Do not use technical jargon unless you explain it.
"""
    response = groq_client.chat.completions.create(
        model="llama-3.3-70b-versatile",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.3
    )
    return response.choices[0].message.content.strip()

@app.route('/voice/query', methods=['POST', 'OPTIONS'])
def voice_query():
    if request.method == 'OPTIONS':
        return make_response()
    data = request.get_json()
    user_query = data.get('query', '').strip()
    sequence = data.get('sequence', '').upper().strip()
    context = data.get('context', '')

    if not user_query:
        return jsonify({"error": "No query provided"}), 400

    try:
        intent_data = parse_intent(user_query, sequence)
    except Exception as e:
        return jsonify({"error": f"Failed to parse intent: {str(e)}"}), 500

    intent = intent_data.get('intent')
    result = {}

    if intent == 'pattern_search' and sequence and intent_data.get('pattern'):
        pattern = intent_data['pattern'].upper().strip()
        positions = kmp_search(sequence, pattern)
        gc = calculate_gc(sequence)
        result = {"pattern": pattern, "occurrences": len(positions), "positions": positions, "gc_content": gc}
    elif intent == 'mutation_detection' and sequence and intent_data.get('second_sequence'):
        reference = intent_data['second_sequence'].upper().strip()
        mutations = detect_mutations(reference, sequence)
        result = {"total_mutations": len(mutations), "mutations": mutations}
    elif intent == 'alignment' and sequence and intent_data.get('second_sequence'):
        seq2 = intent_data['second_sequence'].upper().strip()
        score, aligned1, aligned2, similarity = smith_waterman(sequence, seq2)
        result = {"alignment_score": score, "similarity_percentage": similarity, "aligned_seq1": aligned1, "aligned_seq2": aligned2}
    elif intent == 'gc_content' and sequence:
        gc = calculate_gc(sequence)
        result = {"gc_content": gc}
    elif intent == 'explain':
        result = {"topic": intent_data.get('plain_english_query', user_query)}
    else:
        result = {"note": "Could not determine analysis type from query"}

    try:
        explanation = generate_plain_english_response(intent, result, user_query, context)
    except Exception as e:
        explanation = "Analysis complete. Please see the results below."

    return jsonify({"intent": intent, "result": result, "explanation": explanation})

if __name__ == '__main__':
    app.run(debug=True)