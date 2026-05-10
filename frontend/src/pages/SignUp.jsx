import { useNavigate } from 'react-router-dom'
import { useState } from 'react'

export default function SignUp() {
  const navigate = useNavigate()
  const [form, setForm] = useState({ username: '', email: '', password: '', confirm: '' })

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-900 via-blue-800 to-teal-500 flex flex-col items-center justify-center">
      <div className="bg-white rounded-2xl p-8 w-full max-w-md shadow-2xl">
        <div className="mb-4">
          <label className="text-gray-700 font-medium text-sm mb-1 block">Username</label>
          <div className="flex items-center border border-gray-200 rounded-xl px-4 py-3 gap-2">
            <span className="text-gray-400">👤</span>
            <input className="flex-1 outline-none text-gray-700" placeholder="Choose a username"
              value={form.username} onChange={e => setForm({...form, username: e.target.value})} />
          </div>
        </div>
        <div className="mb-4">
          <label className="text-gray-700 font-medium text-sm mb-1 block">Email Address</label>
          <div className="flex items-center border border-gray-200 rounded-xl px-4 py-3 gap-2">
            <span className="text-gray-400">✉️</span>
            <input className="flex-1 outline-none text-gray-700" placeholder="Enter your email"
              value={form.email} onChange={e => setForm({...form, email: e.target.value})} />
          </div>
        </div>
        <div className="mb-4">
          <label className="text-gray-700 font-medium text-sm mb-1 block">Password</label>
          <div className="flex items-center border border-gray-200 rounded-xl px-4 py-3 gap-2">
            <span className="text-gray-400">🔒</span>
            <input type="password" className="flex-1 outline-none text-gray-700" placeholder="Create a password"
              value={form.password} onChange={e => setForm({...form, password: e.target.value})} />
          </div>
          <p className="text-gray-400 text-xs mt-1">Must be at least 6 characters</p>
        </div>
        <div className="mb-6">
          <label className="text-gray-700 font-medium text-sm mb-1 block">Confirm Password</label>
          <div className="flex items-center border border-gray-200 rounded-xl px-4 py-3 gap-2">
            <span className="text-gray-400">🔒</span>
            <input type="password" className="flex-1 outline-none text-gray-700" placeholder="Confirm your password"
              value={form.confirm} onChange={e => setForm({...form, confirm: e.target.value})} />
          </div>
        </div>
        <button onClick={() => navigate('/dashboard')}
          className="w-full bg-gradient-to-r from-blue-900 to-teal-500 text-white font-bold py-3 rounded-xl hover:opacity-90">
          Create Account
        </button>
        <div className="border-t border-gray-100 mt-4 pt-4 text-center">
          <span className="text-gray-500 text-sm">Already have an account? </span>
          <button onClick={() => navigate('/signin')} className="text-blue-900 font-bold text-sm">Sign In</button>
        </div>
      </div>
      <div className="bg-white/10 rounded-xl px-6 py-3 mt-4 flex items-center gap-2">
        <span className="text-white">✓</span>
        <p className="text-white/80 text-sm">Your data is encrypted and secure. We follow medical-grade security standards.</p>
      </div>
    </div>
  )
}