import { useNavigate } from 'react-router-dom'

export default function Dashboard() {
  const navigate = useNavigate()
  return (
    <div className="min-h-screen bg-white">
      <div className="flex items-center justify-between px-10 py-6 border-b border-gray-100">
        <div>
          <h1 className="text-3xl font-black text-blue-900">VoiceGen</h1>
          <p className="text-gray-400 text-sm">AI-powered DNA analysis and insights</p>
        </div>
      </div>

      <div className="flex items-start justify-between px-10 py-12 max-w-7xl mx-auto gap-16">
        {/* Left — feature cards */}
        <div className="flex-1 flex flex-col gap-4">
          {[
            { icon: "⚡", title: "Fast Analysis", desc: "Get results in seconds with our advanced AI models" },
            { icon: "🛡️", title: "Secure & Private", desc: "Your data is encrypted and never shared" },
            { icon: "🎯", title: "High Accuracy", desc: "Medical-grade analysis with 98%+ confidence" },
          ].map((f, i) => (
            <div key={i} className="flex items-center gap-4 bg-gray-50 rounded-2xl p-5 border border-gray-100">
              <div className="w-12 h-12 bg-teal-50 rounded-xl flex items-center justify-center text-xl flex-shrink-0">{f.icon}</div>
              <div>
                <h3 className="font-bold text-gray-900">{f.title}</h3>
                <p className="text-gray-500 text-sm">{f.desc}</p>
              </div>
            </div>
          ))}

          <button onClick={() => navigate('/upload')}
            className="mt-4 bg-blue-900 text-white font-bold py-4 rounded-xl flex items-center justify-center gap-2 hover:bg-blue-800">
            ↑ Upload DNA Sequence
          </button>
        </div>

        {/* Right — mic button */}
        <div className="flex flex-col items-center justify-center gap-6">
          <button onClick={() => navigate('/upload')}
            className="w-56 h-56 rounded-full bg-gradient-to-br from-blue-900 to-teal-400 flex items-center justify-center shadow-2xl hover:scale-105 transition-transform">
            <span className="text-white text-7xl">🎤</span>
          </button>
          <p className="text-gray-400 text-sm">Click to speak or upload a file</p>
        </div>
      </div>
    </div>
  )
}