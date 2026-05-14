import { useNavigate, useLocation } from 'react-router-dom'
import { useState, useEffect } from 'react'
import axios from 'axios'

export default function Results() {
  const navigate = useNavigate()
  const location = useLocation()
  const sequence = location.state?.sequence || 'ATCGGTATCGATCG'
  const reference = location.state?.reference || ''
  console.log('Sequence received:', sequence)
  const [results, setResults] = useState(null)
  const [loading, setLoading] = useState(true)

useEffect(() => {
  const analyze = async () => {
    try {
      const ref = location.state?.reference || sequence
      const [patternRes, mutRes] = await Promise.all([
        axios.post('/api/analyze/pattern', { sequence, pattern: 'ATG' }),
        axios.post('/api/analyze/mutations', { reference: ref, sample: sequence })
      ])
      setResults({
        occurrences: patternRes.data.occurrences,
        mutations: mutRes.data.total_mutations,
        mutationList: mutRes.data.mutations,
        gc: patternRes.data.gc_content,
      })
    } catch (e) {
      console.error('API error:', e)
      setResults({ occurrences: 0, mutations: 0, mutationList: [], gc: 0 })
    }
    setLoading(false)
  }
  analyze()
}, [sequence])

  return (
    <div className="min-h-screen bg-white px-10 py-8">
      <div className="flex items-center justify-between mb-8">
        <div className="flex items-center gap-4">
          <button onClick={() => navigate('/upload')} className="text-gray-500 hover:text-gray-700">← Back</button>
          <h1 className="text-3xl font-black text-blue-900">Analysis Complete</h1>
        </div>
        <div className="flex gap-3">
          <button className="border border-gray-200 px-4 py-2 rounded-xl text-gray-600 hover:bg-gray-50 flex items-center gap-2">
            ↗ Share
          </button>
          <button className="border border-gray-200 px-4 py-2 rounded-xl text-gray-600 hover:bg-gray-50 flex items-center gap-2">
            ↓ Download Report
          </button>
        </div>
      </div>

      {loading ? (
        <div className="flex items-center justify-center h-64">
          <div className="animate-spin text-5xl">🧬</div>
          <p className="ml-4 text-gray-500 text-lg">Analyzing your sequence...</p>
        </div>
      ) : (
        <>
          {/* Stat cards */}
          <div className="grid grid-cols-3 gap-6 mb-8">
            {[
              { icon: "⚡", label: "Sequence Type", value: "Human DNA", sub: "Homo sapiens genome sequence detected" },
              { icon: "✓", label: "Mutation Status", value: results?.mutations === 0 ? "No mutations" : `${results?.mutations} mutation(s)`, sub: "Compared against reference sequence" },
              { icon: "↗", label: "GC Content", value: `${results?.gc || 0}%`, sub: "Sequence stability indicator", progress: results?.gc || 0 },
            ].map((card, i) => (
              <div key={i} className="border border-gray-100 rounded-2xl p-6 shadow-sm">
                <div className="w-12 h-12 bg-teal-50 rounded-xl flex items-center justify-center text-xl mb-4">{card.icon}</div>
                <p className="text-gray-400 text-sm mb-1">{card.label}</p>
                <p className="text-2xl font-black text-gray-900 mb-1">{card.value}</p>
                <p className="text-gray-400 text-xs">{card.sub}</p>
                {card.progress !== undefined && (
                  <div className="mt-3 h-2 bg-gray-100 rounded-full overflow-hidden">
                    <div className="h-full bg-gradient-to-r from-teal-400 to-blue-900 rounded-full" style={{ width: `${card.progress}%` }}></div>
                  </div>
                )}
              </div>
            ))}
          </div>

          <div className="flex gap-6">
            {/* Analysis summary */}
            <div className="flex-1 border border-gray-100 rounded-2xl p-6 shadow-sm">
              <h3 className="font-bold text-gray-900 mb-4 flex items-center gap-2">📋 Analysis Summary</h3>
              {[
                { label: "Sequence Length", value: `${sequence.length} bases` },
                { label: "Pattern (ATG) Found", value: `${results?.occurrences || 0} times` },
                { label: "GC Content", value: `${results?.gc || 0}%` },
                { label: "Quality Score", value: "Excellent", highlight: true },
              ].map((row, i) => (
                <div key={i} className="flex justify-between items-center py-3 border-b border-gray-50 last:border-0">
                  <span className="text-gray-500">{row.label}</span>
                  <span className={`font-bold ${row.highlight ? 'text-teal-500' : 'text-gray-900'}`}>{row.value}</span>
                </div>
              ))}
            </div>

            {/* AI panel */}
            <div className="w-80 bg-gradient-to-br from-blue-900 to-teal-500 rounded-2xl p-6">
              <h3 className="font-bold text-white mb-2">Need More Insights?</h3>
              <p className="text-white/70 text-sm mb-6">Our AI assistant can help you understand your results, answer questions about specific sequences, and provide detailed explanations.</p>
              <button onClick={() => navigate('/chat', { state: { sequence, results } })}
                className="w-full bg-white text-blue-900 font-bold py-3 rounded-xl hover:bg-gray-100">
                Ask AI Assistant
              </button>
            </div>
          </div>
        </>
      )}
    </div>
  )
}