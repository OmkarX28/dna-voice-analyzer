import { useNavigate } from 'react-router-dom'
import { useState, useRef } from 'react'

const DEMO_SEQUENCES = {
  hbb: {
    label: "🧬 HBB — Sickle Cell",
    condition: "Sickle Cell Anemia",
    sample:    "ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTGTGAAGTCTGCCGTTACTGCC",
    reference: "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTGTGAAGTCTGCCGTTACTGCC",
  },
  brca1: {
    label: "🎗️ BRCA1 — Breast Cancer",
    condition: "BRCA1 Breast Cancer Risk",
    sample:    "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAGGTAATAGAGTCCCGGGAAAGGGACAGGGGGCCCAAGTGATGCTCTG",
    reference: "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGCTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAGGTAATAGAGTCCCGGGAAAGGGACAGGGGGCCCAAGTGATGCTCTG",
  },
  brca2: {
    label: "🎗️ BRCA2 — Breast Cancer",
    condition: "BRCA2 Breast Cancer Risk",
    sample:    "AGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGGTCGCGCCGGGAGAAGCGTGAGGGGACAGATTTGTGACCGGCGCGGTTTTTGTCAGCTTACTCCGGCCAAAAAA",
    reference: "AGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGGTCGCGCCGGGAGAAGCGTGAGGGGACAGATTTGTGACCGGCGCGGTTTTTGTCAGCTTACTCCGGCCAAAAAAT",
  },
  cftr: {
    label: "🫁 CFTR — Cystic Fibrosis",
    condition: "Cystic Fibrosis",
    sample:    "GTAGTAGGTCTTTGGCATTAGGAGCTTGAGCCCAGACGGCCCTAGCAGGGACCCCAGCGCCCGAGAGACCATGCAGAGGTCGCCTCTGGAAAAGGCCAGCGTTGTCTCCAAACTTTTTTTCAGCTGGACCAGACCAATTT",
    reference: "GTAGTAGGTCTTTGGCATTAGGAGCTTGAGCCCAGACGGCCCTAGCAGGGACCCCAGCGCCCGAGAGACCATGCAGAGGTCGCCTCTGGAAAAGGCCAGCGTTTTCTCCAAACTTTTTTTCAGCTGGACCAGACCAATTT",
  },
  jak2: {
    label: "🩸 JAK2 — Blood Cancer",
    condition: "JAK2 Blood Cancer Mutation",
    sample:    "ATTCGGGGAGACTGCAGGCCAACCGGGAGGCTGAGTTCGAAGCTAGCAGGGCGGCGAAGCCAGTGTCGCCCGCGGCGCGTTGAGAAGACGGTGTGGCCCCCGGAGAGGGGTGGAGACAACTGTGACGGGCTTGCCCG",
    reference: "ATTCGGGGAGACTGCAGGCCAACCGGGAGGCTGAGTTCGAAGCTAGCAGGGCGGCGAAGCCAGTGTCGCCCGCGGCGCGTTGAGAAGACGGTGTGGCCCCCGGAGAGGGGTGGAGACAACTGTGACGGGCTTGCCCA",
  },
  kras: {
    label: "🔬 KRAS — Lung Cancer",
    condition: "KRAS Lung/Colorectal Cancer",
    sample:    "CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCAGTGGCGGCGGCGGCGAAGGTGGCGGGCTCGGCCAGTACTCCCGGCCCCGCCATTTCGGACTGGGAGCGAGCGCGGCGCAGGCACTGAAGGCGGCGGC",
    reference: "CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCAGTGGCGGCGGCGGCGAAGGTGGCGGGCTCGGCCAGTACTCCCGGCCCCGCCATTTCGGACTGGGAGCGAGCGCGGCGCAGGCACTGAAGGCGGCGGT",
  },
}

export default function Upload() {
  const navigate = useNavigate()
  const [sequence, setSequence] = useState('')
  const [reference, setReference] = useState('')
  const [activeDemo, setActiveDemo] = useState(null)
  const [dragging, setDragging] = useState(false)
  const fileRef = useRef()

  const handleFile = (file) => {
    const reader = new FileReader()
    reader.onload = (e) => {
      let text = e.target.result
      const lines = text.split('\n').filter(l => !l.startsWith('>')).join('')
      setSequence(lines.toUpperCase().replace(/[^ATCG]/g, ''))
    }
    reader.readAsText(file)
  }

  const loadDemo = (key) => {
    const demo = DEMO_SEQUENCES[key]
    setSequence(demo.sample)
    setReference(demo.reference)
    setActiveDemo(key)
  }

  return (
    <div className="min-h-screen bg-white px-10 py-8">
      <button onClick={() => navigate('/dashboard')} className="flex items-center gap-2 text-gray-500 hover:text-gray-700 mb-6">
        ← Back
      </button>
      <h1 className="text-3xl font-black text-blue-900 mb-8">Upload DNA Sequence</h1>

      <div className="flex gap-8">
        <div className="flex-1">
          {/* Drag and drop */}
          <div
            onDragOver={e => { e.preventDefault(); setDragging(true) }}
            onDragLeave={() => setDragging(false)}
            onDrop={e => { e.preventDefault(); setDragging(false); handleFile(e.dataTransfer.files[0]) }}
            className={`border-2 border-dashed rounded-2xl p-12 flex flex-col items-center justify-center gap-4 transition-colors ${dragging ? 'border-teal-400 bg-teal-50' : 'border-gray-200 bg-gray-50'}`}>
            <div className="w-16 h-16 bg-blue-100 rounded-full flex items-center justify-center">
              <span className="text-blue-900 text-3xl">↑</span>
            </div>
            <p className="text-gray-700 font-semibold text-lg">Drag and drop your FASTA file here</p>
            <p className="text-gray-400">or</p>
            <button onClick={() => fileRef.current.click()}
              className="bg-blue-900 text-white font-bold px-8 py-3 rounded-xl hover:bg-blue-800">
              Browse Files
            </button>
            <input ref={fileRef} type="file" accept=".txt,.fasta,.fa,.seq,.fna" className="hidden"
              onChange={e => handleFile(e.target.files[0])} />
            <p className="text-gray-400 text-sm">Supported formats: .txt, .fasta, .fa, .seq, .fna</p>
          </div>

          {/* Sample sequence */}
          <div className="mt-6">
            <label className="text-gray-700 font-semibold mb-2 block">
              Sample Sequence <span className="text-red-400">*</span>
              <span className="text-gray-400 font-normal text-sm ml-2">— the sequence you want to analyze</span>
            </label>
            <textarea
              className="w-full border border-gray-200 rounded-xl p-4 text-gray-700 font-mono text-sm h-28 outline-none focus:border-teal-400"
              placeholder="Paste your DNA sequence here..."
              value={sequence}
              onChange={e => { setSequence(e.target.value.toUpperCase().replace(/[^ATCG\n]/g, '')); setActiveDemo(null) }}
            />
          </div>

          {/* Reference sequence */}
          <div className="mt-4">
            <label className="text-gray-700 font-semibold mb-2 block">
              Reference Sequence
              <span className="text-gray-400 font-normal text-sm ml-2">— optional, used for mutation detection</span>
            </label>
            <textarea
              className="w-full border border-gray-200 rounded-xl p-4 text-gray-700 font-mono text-sm h-28 outline-none focus:border-teal-400"
              placeholder="Paste a reference sequence to compare against..."
              value={reference}
              onChange={e => { setReference(e.target.value.toUpperCase().replace(/[^ATCG\n]/g, '')); setActiveDemo(null) }}
            />
          </div>

          {/* Real gene demo buttons */}
          <div className="mt-4 bg-blue-50 rounded-xl p-4">
            <p className="text-blue-900 font-semibold text-sm mb-3">⚡ Quick Load — Real Gene Sequences from NCBI:</p>
            <div className="grid grid-cols-2 gap-2">
              {Object.entries(DEMO_SEQUENCES).map(([key, demo]) => (
                <button
                  key={key}
                  onClick={() => loadDemo(key)}
                  className={`text-left px-4 py-3 rounded-xl text-sm font-medium border transition-all ${
                    activeDemo === key
                      ? 'bg-blue-900 text-white border-blue-900'
                      : 'bg-white border-blue-200 text-blue-900 hover:bg-blue-100'
                  }`}>
                  {demo.label}
                  <p className={`text-xs mt-0.5 font-normal ${activeDemo === key ? 'text-blue-200' : 'text-gray-400'}`}>
                    {demo.condition}
                  </p>
                </button>
              ))}
            </div>
          </div>

          <button
            onClick={() => navigate('/results', { state: { sequence, reference, condition: activeDemo ? DEMO_SEQUENCES[activeDemo].condition : null } })}
            disabled={!sequence}
            className="mt-6 w-full bg-gradient-to-r from-blue-900 to-teal-500 text-white font-bold py-4 rounded-xl hover:opacity-90 disabled:opacity-40 disabled:cursor-not-allowed">
            Analyze Sequence →
          </button>
        </div>

        {/* Right panel */}
        <div className="w-72 flex flex-col gap-4">
          <div className="bg-white border border-gray-100 rounded-2xl p-6 shadow-sm">
            <h3 className="font-bold text-gray-900 mb-4">Upload Guidelines</h3>
            {[
              { icon: "📄", text: "Supports FASTA, .fna, plain text formats" },
              { icon: "⏱️", text: "Analysis completes within 10-30 seconds" },
              { icon: "🧬", text: "Provide a reference sequence for mutation detection" },
              { icon: "⚡", text: "Use Quick Load for real NCBI gene sequences" },
            ].map((g, i) => (
              <div key={i} className="flex items-start gap-3 mb-3">
                <span className="text-teal-500">{g.icon}</span>
                <p className="text-gray-500 text-sm">{g.text}</p>
              </div>
            ))}
          </div>

          <div className="bg-gradient-to-br from-blue-900 to-teal-500 rounded-2xl p-6">
            <h3 className="font-bold text-white mb-2">Gene Database</h3>
            <p className="text-white/70 text-sm">All sequences sourced from NCBI — America's National Center for Biotechnology Information. Medically validated by our Biotech team.</p>
          </div>
        </div>
      </div>
    </div>
  )
}