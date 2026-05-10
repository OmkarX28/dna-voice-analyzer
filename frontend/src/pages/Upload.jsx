import { useNavigate } from 'react-router-dom'
import { useState, useRef } from 'react'

export default function Upload() {
  const navigate = useNavigate()
  const [sequence, setSequence] = useState('')
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
          <div
            onDragOver={e => { e.preventDefault(); setDragging(true) }}
            onDragLeave={() => setDragging(false)}
            onDrop={e => { e.preventDefault(); setDragging(false); handleFile(e.dataTransfer.files[0]) }}
            className={`border-2 border-dashed rounded-2xl p-16 flex flex-col items-center justify-center gap-4 transition-colors ${dragging ? 'border-teal-400 bg-teal-50' : 'border-gray-200 bg-gray-50'}`}>
            <div className="w-16 h-16 bg-blue-100 rounded-full flex items-center justify-center">
              <span className="text-blue-900 text-3xl">↑</span>
            </div>
            <p className="text-gray-700 font-semibold text-lg">Drag and drop your file here</p>
            <p className="text-gray-400">or</p>
            <button onClick={() => fileRef.current.click()}
              className="bg-blue-900 text-white font-bold px-8 py-3 rounded-xl hover:bg-blue-800">
              Browse Files
            </button>
            <input ref={fileRef} type="file" accept=".txt,.fasta,.fa,.seq" className="hidden"
              onChange={e => handleFile(e.target.files[0])} />
            <p className="text-gray-400 text-sm">Supported formats: .txt, .fasta, .fa, .seq</p>
          </div>

          {/* Manual input */}
          <div className="mt-6">
            <label className="text-gray-700 font-semibold mb-2 block">Or paste sequence directly:</label>
            <textarea
              className="w-full border border-gray-200 rounded-xl p-4 text-gray-700 font-mono text-sm h-32 outline-none focus:border-teal-400"
              placeholder="Paste your DNA sequence here (e.g. ATCGGTATCGATCG...)"
              value={sequence}
              onChange={e => setSequence(e.target.value.toUpperCase().replace(/[^ATCG\n]/g, ''))}
            />
          </div>

          <button
            onClick={() => navigate('/results', { state: { sequence } })}
            disabled={!sequence}
            className="mt-4 w-full bg-gradient-to-r from-blue-900 to-teal-500 text-white font-bold py-4 rounded-xl hover:opacity-90 disabled:opacity-40 disabled:cursor-not-allowed">
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
              { icon: "↑", text: "Maximum file size: 50MB" },
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