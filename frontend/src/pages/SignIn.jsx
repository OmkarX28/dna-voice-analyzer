import { useNavigate } from 'react-router-dom'
import { useState } from 'react'

export default function SignIn() {
  const navigate = useNavigate()
  const [form, setForm] = useState({ email: '', password: '' })

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-900 via-blue-800 to-teal-500 flex flex-col items-center justify-center">
      <div className="w-16 h-16 rounded-full bg-blue-800 flex items-center justify-center mb-4">
        <span className="text-white text-2xl">🧬</span>
      </div>
      <h1 className="text-white text-3xl font-black mb-1">VoiceGen</h1>
      <p className="text-white/70 mb-8">Sign in to access your DNA analysis platform</p>

      <div className="bg-white rounded-2xl p-8 w-full max-w-md shadow-2xl">
        <div className="mb-4">
          <label className="text-gray-700 font-medium text-sm mb-1 block">Username or Email</label>
          <div className="flex items-center border border-gray-200 rounded-xl px-4 py-3 gap-2">
            <span className="text-gray-400">✉️</span>
            <input className="flex-1 outline-none text-gray-700" placeholder="Enter username or email"
              value={form.email} onChange={e => setForm({...form, email: e.target.value})} />
          </div>
        </div>
        <div className="mb-6">
          <label className="text-gray-700 font-medium text-sm mb-1 block">Password</label>
          <div className="flex items-center border border-gray-200 rounded-xl px-4 py-3 gap-2">
            <span className="text-gray-400">🔒</span>
            <input type="password" className="flex-1 outline-none text-gray-700" placeholder="Enter password"
              value={form.password} onChange={e => setForm({...form, password: e.target.value})} />
          </div>
        </div>
        <button onClick={() => navigate('/dashboard')}
          className="w-full bg-gradient-to-r from-blue-900 to-teal-500 text-white font-bold py-3 rounded-xl hover:opacity-90">
          Sign In
        </button>
        <p className="text-center text-gray-500 text-sm mt-4">
          Demo: admin@voicegen.com / dna123456
        </p>
        <div className="border-t border-gray-100 mt-4 pt-4 text-center">
          <span className="text-gray-500 text-sm">Don't have an account? </span>
          <button onClick={() => navigate('/signup')} className="text-blue-900 font-bold text-sm">Create Account</button>
        </div>
      </div>
      <p className="text-white/50 text-sm mt-6">Secure medical-grade authentication powered by AI</p>
    </div>
  )
}