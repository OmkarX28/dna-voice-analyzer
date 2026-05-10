import { useNavigate } from 'react-router-dom'

export default function Landing() {
  const navigate = useNavigate()
  return (
    <div className="min-h-screen bg-white">
      {/* Navbar */}
      <nav className="flex items-center justify-between px-10 py-4 border-b border-gray-100">
        <div className="flex items-center gap-2">
          <div className="w-9 h-9 rounded-lg bg-gradient-to-br from-blue-900 to-teal-500 flex items-center justify-center">
            <span className="text-white text-sm font-bold">V</span>
          </div>
          <span className="text-xl font-bold text-gray-900">VoiceGen</span>
        </div>
        <div className="flex items-center gap-4">
          <button onClick={() => navigate('/signin')} className="text-gray-600 font-medium hover:text-blue-900">Sign In</button>
          <button onClick={() => navigate('/signup')} className="bg-blue-900 text-white px-5 py-2 rounded-lg font-medium hover:bg-blue-800">Get Started</button>
        </div>
      </nav>

      {/* Hero */}
      <div className="flex items-center justify-between px-10 py-20 max-w-7xl mx-auto">
        <div className="max-w-xl">
          <span className="bg-teal-50 text-teal-600 text-sm px-4 py-1 rounded-full font-medium">AI-Powered DNA Analysis</span>
          <h1 className="text-5xl font-black text-gray-900 mt-4 leading-tight">
            Transform DNA<br />Analysis with{' '}
            <span className="text-blue-900">Voice</span>{' '}
            <span className="text-teal-500">& AI</span>
          </h1>
          <p className="text-gray-500 mt-4 text-lg leading-relaxed">
            VoiceGen combines cutting-edge voice recognition with advanced AI to deliver fast, accurate DNA sequence analysis. Simply speak or upload your data and get medical-grade insights in seconds.
          </p>
          <div className="flex gap-4 mt-8">
            <button onClick={() => navigate('/signup')} className="bg-blue-900 text-white px-6 py-3 rounded-lg font-semibold hover:bg-blue-800 flex items-center gap-2">
              Get Started Free →
            </button>
            <button onClick={() => navigate('/signin')} className="border border-gray-300 text-gray-700 px-6 py-3 rounded-lg font-semibold hover:bg-gray-50">
              Sign In
            </button>
          </div>
          <div className="flex gap-10 mt-10">
            <div><p className="text-2xl font-black text-blue-900">98.7%</p><p className="text-gray-400 text-sm">Accuracy Rate</p></div>
            <div className="border-l border-gray-200 pl-10"><p className="text-2xl font-black text-blue-900">10k+</p><p className="text-gray-400 text-sm">Analyses Complete</p></div>
            <div className="border-l border-gray-200 pl-10"><p className="text-2xl font-black text-blue-900">&lt;30s</p><p className="text-gray-400 text-sm">Average Time</p></div>
          </div>
        </div>

        {/* Hero visual */}
        <div className="w-96 h-72 bg-gradient-to-br from-blue-900 to-teal-400 rounded-2xl p-6 flex flex-col justify-between shadow-2xl">
          <div className="bg-white/20 rounded-xl p-4 flex items-center gap-3">
            <div className="w-10 h-10 bg-white/30 rounded-full flex items-center justify-center">
              <span className="text-white text-lg">🎤</span>
            </div>
            <div className="flex-1">
              <div className="h-2 bg-white/40 rounded w-3/4 mb-2"></div>
              <div className="h-2 bg-white/30 rounded w-1/2"></div>
            </div>
          </div>
          <div className="flex gap-3">
            <div className="flex-1 bg-white/20 rounded-xl p-3">
              <span className="text-white text-lg">✓</span>
              <p className="text-white text-sm font-semibold mt-1">No Mutations</p>
            </div>
            <div className="flex-1 bg-white/20 rounded-xl p-3">
              <span className="text-white text-lg">↗</span>
              <p className="text-white text-sm font-semibold mt-1">98.7% Confidence</p>
            </div>
          </div>
        </div>
      </div>

      {/* Features Grid */}
      <div className="bg-gray-50 px-10 py-16">
        <div className="max-w-7xl mx-auto grid grid-cols-3 gap-6">
          {[
            { icon: "🎤", title: "Voice-Powered Analysis", desc: "Speak your DNA sequences naturally. Our advanced voice recognition understands scientific terminology." },
            { icon: "⚡", title: "Real-Time Processing", desc: "Get results in under 30 seconds. Our AI analyzes millions of base pairs instantly with medical-grade accuracy." },
            { icon: "🛡️", title: "Bank-Level Security", desc: "Your data is encrypted end-to-end. We comply with HIPAA and follow strict medical data protection standards." },
            { icon: "🤖", title: "AI-Powered Insights", desc: "Chat with our AI assistant to understand complex results. Get explanations in plain language or technical detail." },
            { icon: "🧬", title: "Mutation Detection", desc: "Automatically identify genetic mutations and variations. Receive detailed reports with confidence scores." },
            { icon: "📊", title: "Detailed Analytics", desc: "Access comprehensive reports with visualizations, quality scores, and downloadable results for your records." },
          ].map((f, i) => (
            <div key={i} className="bg-white rounded-2xl p-6 shadow-sm border border-gray-100">
              <div className="w-12 h-12 bg-teal-50 rounded-xl flex items-center justify-center text-2xl mb-4">{f.icon}</div>
              <h3 className="font-bold text-gray-900 mb-2">{f.title}</h3>
              <p className="text-gray-500 text-sm leading-relaxed">{f.desc}</p>
            </div>
          ))}
        </div>
      </div>

      {/* CTA */}
      <div className="bg-gradient-to-r from-blue-900 to-teal-500 py-20 text-center">
        <h2 className="text-4xl font-black text-white mb-4">Ready to Transform Your DNA Analysis?</h2>
        <p className="text-white/80 mb-8">Join thousands of healthcare professionals using VoiceGen for faster, more accurate genetic analysis.</p>
        <button onClick={() => navigate('/signup')} className="bg-white text-blue-900 font-bold px-8 py-3 rounded-xl hover:bg-gray-100">
          Get Started Free
        </button>
      </div>

      {/* Footer */}
      <div className="bg-gray-900 px-10 py-6 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-blue-900 to-teal-500 flex items-center justify-center">
            <span className="text-white text-xs font-bold">V</span>
          </div>
          <span className="text-white font-bold">VoiceGen</span>
        </div>
        <p className="text-gray-400 text-sm">© 2026 VoiceGen. All rights reserved.</p>
      </div>
    </div>
  )
}