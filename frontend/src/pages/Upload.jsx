import { useNavigate } from 'react-router-dom'
import { useState, useRef } from 'react'

export default function Upload() {
  const navigate = useNavigate()
  const [sequence, setSequence] = useState('')
  const [reference, setReference] = useState('')
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

  return (
    <div className="min-h-screen bg-white px-10 py-8">
      <button onClick={() => navigate('/dashboard')} className="flex items-center gap-2 text-gray-500 hover:text-gray-700 mb-6">
        ← Back
      </button>
      <h1 className="text-3xl font-black text-blue-900 mb-8">Upload DNA Sequence</h1>

      <div className="flex gap-8">
        {/* Upload area */}
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
            <input ref={fileRef} type="file" accept=".txt,.fasta,.fa,.seq" className="hidden"
              onChange={e => handleFile(e.target.files[0])} />
            <p className="text-gray-400 text-sm">Supported formats: .txt, .fasta, .fa, .seq</p>
          </div>

          {/* Sample sequence input */}
          <div className="mt-6">
            <label className="text-gray-700 font-semibold mb-2 block">
              Sample Sequence <span className="text-red-400">*</span>
              <span className="text-gray-400 font-normal text-sm ml-2">— the sequence you want to analyze</span>
            </label>
            <textarea
              className="w-full border border-gray-200 rounded-xl p-4 text-gray-700 font-mono text-sm h-28 outline-none focus:border-teal-400"
              placeholder="Paste your DNA sequence here (e.g. ATGGTGCACCTGACT...)"
              value={sequence}
              onChange={e => setSequence(e.target.value.toUpperCase().replace(/[^ATCG\n]/g, ''))}
            />
          </div>

          {/* Reference sequence input */}
          <div className="mt-4">
            <label className="text-gray-700 font-semibold mb-2 block">
              Reference Sequence
              <span className="text-gray-400 font-normal text-sm ml-2">— optional, used for mutation detection</span>
            </label>
            <textarea
              className="w-full border border-gray-200 rounded-xl p-4 text-gray-700 font-mono text-sm h-28 outline-none focus:border-teal-400"
              placeholder="Paste a reference sequence to compare against (e.g. normal BRCA1 gene snippet)"
              value={reference}
              onChange={e => setReference(e.target.value.toUpperCase().replace(/[^ATCG\n]/g, ''))}
            />
          </div>

          {/* Quick demo sequences */}
          <div className="mt-4 bg-blue-50 rounded-xl p-4">
            <p className="text-blue-900 font-semibold text-sm mb-3">Quick Demo — load a real mutation example:</p>
            <div className="flex gap-3 flex-wrap">
              <button
                onClick={() => {
                  setSequence('ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG')
                  setReference('ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG')
                }}
                className="bg-white border border-blue-200 text-blue-900 text-xs px-4 py-2 rounded-lg hover:bg-blue-100 font-medium">
                🧬 Sickle Cell Mutation
              </button>
              <button
                onClick={() => {
                  setSequence('ATCGGTATCGATCGGCTAGCTAGCAAAC')
                  setReference('ATCGGTATCGATCGGCTAGCTAGCTAGC')
                }}
                className="bg-white border border-blue-200 text-blue-900 text-xs px-4 py-2 rounded-lg hover:bg-blue-100 font-medium">
                🔬 Multiple Point Mutations
              </button>
              <button
                onClick={() => {
                  setSequence('ATCGGTATCGATCGTTTT')
                  setReference('ATCGGTATCGATCG')
                }}
                className="bg-white border border-blue-200 text-blue-900 text-xs px-4 py-2 rounded-lg hover:bg-blue-100 font-medium">
                ➕ Insertion Example
              </button>
            </div>
          </div>

          <button
            onClick={() => navigate('/results', { state: { sequence, reference } })}
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
              { icon: "📄", text: "Ensure your DNA sequence is in FASTA or plain text format" },
              { icon: "⏱️", text: "Analysis typically completes within 10-30 seconds" },
              { icon: "🧬", text: "Provide a reference sequence for mutation detection" },
            ].map((g, i) => (
              <div key={i} className="flex items-start gap-3 mb-3">
                <span className="text-teal-500">{g.icon}</span>
                <p className="text-gray-500 text-sm">{g.text}</p>
              </div>
            ))}
          </div>

          <div className="bg-gradient-to-br from-blue-900 to-teal-500 rounded-2xl p-6">
            <h3 className="font-bold text-white mb-2">Recent Uploads</h3>
            <p className="text-white/60 text-sm">No recent uploads found</p>
            <p className="text-white/40 text-xs mt-2">Your upload history will appear here for quick access</p>
          </div>
        </div>
      </div>
    </div>
  )
}