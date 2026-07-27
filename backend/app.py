from flask import Flask, request, jsonify, make_response
from flask_cors import CORS
from dotenv import load_dotenv
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

# ── Mutation Knowledge Base ─────────────────────────────────
MUTATION_KNOWLEDGE_BASE = """
KNOWN MUTATIONS REFERENCE — use ONLY when the exact base change matches:

1. Sickle Cell Anemia - HBB gene
   TRIGGER: Point mutation GAG→GTG (G changed to T at the codon 6 position, specifically position where reference=G and sample=T in the HBB gene)
   Explanation: The single base change from GAG to GTG causes hemoglobin to form rigid rods inside red blood cells, making them sickle-shaped instead of round. This blocks blood vessels causing severe pain, anemia, and organ damage. It is one of the most common genetic disorders in India.

2. BRCA1 Breast Cancer - BRCA1 gene
   TRIGGER: Any deletion, insertion, or point mutation in BRCA1 sequence
   Explanation: Mutations in BRCA1 disable a tumor suppressor protein that repairs damaged DNA. Without this repair mechanism, cells can become cancerous. Carriers have up to 72% lifetime risk of breast cancer and elevated ovarian cancer risk.

3. BRCA2 Breast Cancer - BRCA2 gene
   TRIGGER: Any deletion, insertion, or point mutation in BRCA2 sequence
   Explanation: Like BRCA1, BRCA2 mutations disable DNA repair. Carriers face up to 69% lifetime breast cancer risk plus elevated prostate and pancreatic cancer risk.

4. Cystic Fibrosis - CFTR gene
   TRIGGER: Any deletion or point mutation in CFTR sequence
   Explanation: Mutations in CFTR disrupt chloride transport causing thick sticky mucus to accumulate in the lungs and digestive system. This leads to progressive lung damage and digestive problems requiring lifelong treatment.

5. JAK2 Blood Cancer - JAK2 gene
   TRIGGER: Any point mutation in JAK2 sequence
   Explanation: JAK2 mutations cause uncontrolled blood cell production leading to myeloproliferative disorders including polycythemia vera, essential thrombocythemia, and myelofibrosis — serious blood cancers requiring ongoing treatment.

6. KRAS Lung/Colorectal Cancer - KRAS gene
   TRIGGER: Any point mutation in KRAS sequence
   Explanation: KRAS mutations are among the most common oncogenic mutations in human cancers, found in approximately 30% of all cancers. They lock the KRAS protein in an active state, driving uncontrolled cell division in lung, colorectal, and pancreatic cancers.

7. Beta Thalassemia - HBB gene
   TRIGGER: Deletion in HBB sequence
   Explanation: Deletions in HBB reduce beta-globin production causing severe anemia requiring regular blood transfusions. Common in Mediterranean, Middle Eastern, and South Asian populations.

8. Huntington's Disease - HTT gene
   TRIGGER: Insertion/repeat expansion in HTT sequence
   Explanation: CAG repeat expansions in HTT create an abnormally long protein that progressively destroys brain cells causing movement disorders, cognitive decline, and psychiatric symptoms typically appearing between ages 30-50.
"""

# ── FASTA Parser ────────────────────────────────────────────
def parse_fasta_input(raw_text):
    """
    Robustly parses raw input that may be:
    1. A plain DNA string (no headers)
    2. A single FASTA entry (one > header)
    3. A multi-FASTA string with reference and sample separated by > headers
    Returns: (sequence, reference) tuple — both clean ATCG strings
    """
    raw_text = raw_text.strip()
    
    if '>' not in raw_text:
        # Plain DNA string — no headers
        clean = ''.join(c for c in raw_text.upper() if c in 'ATCG')
        return clean, None
    
    # Split on > headers
    entries = []
    current_header = None
    current_seq = []
    
    for line in raw_text.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                seq = ''.join(c for c in ''.join(current_seq).upper() if c in 'ATCG')
                entries.append((current_header, seq))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)
    
    # Don't forget the last entry
    if current_seq:
        seq = ''.join(c for c in ''.join(current_seq).upper() if c in 'ATCG')
        entries.append((current_header, seq))
    
    if len(entries) == 0:
        return '', None
    elif len(entries) == 1:
        return entries[0][1], None
    else:
        # Multiple entries — first is reference, second is sample
        # Or detect by header keywords
        ref_seq = None
        sample_seq = None
        
        for header, seq in entries:
            header_lower = header.lower() if header else ''
            if any(k in header_lower for k in ['ref', 'reference', 'normal', 'wild', 'wt']):
                ref_seq = seq
            elif any(k in header_lower for k in ['sample', 'mutant', 'patient', 'mut', 'test']):
                sample_seq = seq
        
        # If no keyword match, use order: first=reference, second=sample
        if ref_seq is None and sample_seq is None:
            ref_seq = entries[0][1]
            sample_seq = entries[1][1]
        elif ref_seq is None:
            ref_seq = entries[0][1]
        elif sample_seq is None:
            sample_seq = entries[1][1]
        
        return sample_seq, ref_seq

# ── KMP Algorithm ───────────────────────────────────────────
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
    if m == 0 or n == 0:
        return positions
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

# ── GC Content ──────────────────────────────────────────────
def calculate_gc(sequence):
    if not sequence:
        return 0
    g = sequence.count('G')
    c = sequence.count('C')
    return round((g + c) / len(sequence) * 100, 2)

# ── Mutation Detector ───────────────────────────────────────
def detect_mutations(reference, sample):
    mutations = []
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
    if len(sample) > len(reference):
        inserted = sample[len(reference):]
        mutations.append({
            "type": "insertion",
            "position": len(reference),
            "inserted_bases": inserted,
            "description": f"Insertion of '{inserted}' after position {len(reference) - 1}"
        })
    if len(sample) < len(reference):
        deleted = reference[len(sample):]
        mutations.append({
            "type": "deletion",
            "position": len(sample),
            "deleted_bases": deleted,
            "description": f"Deletion of '{deleted}' starting at position {len(sample)}"
        })
    return mutations

# ── Smith-Waterman ──────────────────────────────────────────
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

# ── NLP Intent Engine ───────────────────────────────────────
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

# ── FIXED: Generate Plain English Response ──────────────────
def generate_plain_english_response(intent, analysis_result, user_query, context=''):
    """
    CRITICAL FIX: The LLM is forced to ground its answer STRICTLY on the
    algorithm's JSON output. It must NOT default to Sickle Cell or any
    disease unless the exact base change in analysis_result matches.
    """

    # Build a precise mutation summary from the actual algorithm output
    mutation_summary = ""
    if analysis_result.get("mutations") and len(analysis_result["mutations"]) > 0:
        mutation_lines = []
        for m in analysis_result["mutations"]:
            if m["type"] == "point_mutation":
                mutation_lines.append(
                    f"Point mutation at position {m['position']}: "
                    f"reference base '{m['reference_base']}' changed to '{m['sample_base']}'"
                )
            elif m["type"] == "insertion":
                mutation_lines.append(
                    f"Insertion of bases '{m['inserted_bases']}' at position {m['position']}"
                )
            elif m["type"] == "deletion":
                mutation_lines.append(
                    f"Deletion of bases '{m['deleted_bases']}' starting at position {m['position']}"
                )
        mutation_summary = "\n".join(mutation_lines)
    elif analysis_result.get("total_mutations", 0) == 0:
        mutation_summary = "No mutations detected. Sequences are identical."
    
    prompt = f"""
You are a medical AI assistant explaining DNA analysis results.

STRICT RULES — you MUST follow these or your answer is wrong:
1. Base your ENTIRE explanation on the ACTUAL MUTATION DATA provided below. Do not guess or use prior knowledge about what disease this "might" be.
2. Only mention Sickle Cell Anemia if you see EXACTLY: reference_base='G' changed to sample_base='T' (GAG to GTG mutation).
3. Only mention BRCA1/BRCA2 if the context explicitly says BRCA gene.
4. Only mention CFTR/Cystic Fibrosis if the context explicitly says CFTR gene.
5. Only mention JAK2 if the context explicitly says JAK2 gene.
6. If you cannot match the mutation to a known disease from the context, simply describe the exact base change found and say it requires clinical evaluation.
7. NEVER hallucinate a disease name that is not supported by the actual mutation data.

ACTUAL ALGORITHM OUTPUT (ground truth — use ONLY this):
{json.dumps(analysis_result, indent=2)}

EXACT MUTATIONS FOUND BY ALGORITHM:
{mutation_summary}

ADDITIONAL CONTEXT: {context}

USER QUESTION: "{user_query}"

KNOWN MUTATIONS REFERENCE (only use if the base change matches exactly):
{MUTATION_KNOWLEDGE_BASE}

Write a clear, friendly 2-3 sentence explanation grounded strictly in the mutation data above.
State the exact base change found. Then if it matches a known disease, name it and explain it briefly.
If it does not match any known disease, say the mutation requires clinical evaluation.
Do not mention JSON, code, or algorithms.
"""
    response = groq_client.chat.completions.create(
        model="llama-3.3-70b-versatile",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.1
    )
    return response.choices[0].message.content.strip()

# ── Routes ──────────────────────────────────────────────────
@app.route('/ping', methods=['GET'])
def ping():
    return jsonify({"message": "DNA Analyzer backend is running!"})

@app.route('/analyze/pattern', methods=['POST', 'OPTIONS'])
def find_pattern():
    if request.method == 'OPTIONS':
        return make_response()
    data = request.get_json()
    raw_sequence = data.get('sequence', '')
    pattern = data.get('pattern', '').upper().strip()

    sequence, _ = parse_fasta_input(raw_sequence)

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
        "message": f"Pattern '{pattern}' found {len(positions)} time(s)." if positions else f"Pattern '{pattern}' was not found."
    })

@app.route('/analyze/mutations', methods=['POST', 'OPTIONS'])
def find_mutations():
    if request.method == 'OPTIONS':
        return make_response()
    data = request.get_json()
    raw_reference = data.get('reference', '')
    raw_sample = data.get('sample', '')

    # Handle multi-FASTA input — if reference is empty but sample has > headers
    if not raw_reference and '>' in raw_sample:
        sample, reference = parse_fasta_input(raw_sample)
    else:
        sample, _ = parse_fasta_input(raw_sample)
        reference, _ = parse_fasta_input(raw_reference)

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
        "gc_content": {"reference": gc_ref, "sample": gc_sample},
        "message": f"{len(mutations)} mutation(s) detected." if mutations else "No mutations detected. Sequences are identical."
    })

@app.route('/analyze/align', methods=['POST', 'OPTIONS'])
def align_sequences():
    if request.method == 'OPTIONS':
        return make_response()
    data = request.get_json()
    raw_seq1 = data.get('seq1', '')
    raw_seq2 = data.get('seq2', '')

    seq1, _ = parse_fasta_input(raw_seq1)
    seq2, _ = parse_fasta_input(raw_seq2)

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

@app.route('/voice/query', methods=['POST', 'OPTIONS'])
def voice_query():
    if request.method == 'OPTIONS':
        return make_response()
    data = request.get_json()
    user_query = data.get('query', '').strip()
    raw_sequence = data.get('sequence', '')
    context = data.get('context', '')

    sequence, parsed_reference = parse_fasta_input(raw_sequence)

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
        result = {
            "pattern": pattern,
            "occurrences": len(positions),
            "positions": positions,
            "gc_content": gc
        }

    elif intent == 'mutation_detection' and sequence:
        reference = intent_data.get('second_sequence', '')
        if reference:
            reference = reference.upper().strip()
        elif parsed_reference:
            reference = parsed_reference
        if reference:
            mutations = detect_mutations(reference, sequence)
            result = {
                "total_mutations": len(mutations),
                "mutations": mutations,
                "gc_content": calculate_gc(sequence)
            }
        else:
            result = {"note": "No reference sequence provided for mutation detection"}

    elif intent == 'alignment' and sequence and intent_data.get('second_sequence'):
        seq2 = intent_data['second_sequence'].upper().strip()
        score, aligned1, aligned2, similarity = smith_waterman(sequence, seq2)
        result = {
            "alignment_score": score,
            "similarity_percentage": similarity,
            "aligned_seq1": aligned1,
            "aligned_seq2": aligned2
        }

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