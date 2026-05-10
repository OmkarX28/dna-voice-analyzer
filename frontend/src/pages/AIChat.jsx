import { useNavigate, useLocation } from 'react-router-dom'
import { useState, useRef, useEffect } from 'react'
import axios from 'axios'

export default function AIChat() {
  const navigate = useNavigate()
  const location = useLocation()
  const sequence = location.state?.sequence || 'ATCGGTATCGATCG'
  const [messages, setMessages] = useState([
    { role: 'ai', text: "Hello! I've analyzed your DNA sequence. What would you like to know?" }
  ])
  const [input, setInput] = useState('')
  const [loading, setLoading] = useState(false)
  const [listening, setListening] = useState(false)
  const bottomRef = useRef()

  useEffect(() => { bottomRef.current?.scrollIntoView({ behavior: 'smooth' }) }, [messages])

  const speak = (text) => {
    const utterance = new SpeechSynthesisUtterance(text)
    utterance.rate = 0.9
    window.speechSynthesis.speak(utterance)
  }

  const startListening = () => {
    const SpeechRecognition = window.SpeechRecognition || window.webkitSpeechRecognition
    if (!SpeechRecognition) { alert('Voice not supported in this browser. Use Chrome.'); return }
    const recognition = new SpeechRecognition()
    recognition.lang = 'en-US'
    recognition.onstart = () => setListening(true)
    recognition.onend = () => setListening(false)
    recognition.onresult = (e) => {
      const transcript = e.results[0][0].transcript
      setInput(transcript)
    }
    recognition.start()
  }

  const sendMessage = async (text) => {
    const query = text || input
    if (!query.trim()) return
    setMessages(prev => [...prev, { role: 'user', text: query }])
    setInput('')
    setLoading(true)
    try {
      const res = await axios.post('http://localhost:5000/voice/query', { query, sequence })
      const reply = res.data.explanation
      setMessages(prev => [...prev, { role: 'ai', text: reply }])
      speak(reply)
    } catch (e) {
      setMessages(prev => [...prev, { role: 'ai', text: 'Sorry, I encountered an error. Please try again.' }])
    }
    setLoading(false)
  }

  const suggestions = [
    "What does the mutation status indicate?",
    "How accurate is this analysis?",
    "Can you explain the sequence type?",
    "What should I do next?",
  ]

  return (
    <div className="min-h-screen bg-white flex">
      {/* Main chat */}
      <div className="flex-1 flex flex-col">
        {/* Header */}
        <div className="flex items-center justify-between px-8 py-4 border-b border-gray-100">
          <button onClick={() => navigate('/results')} className="text-gray-500 hover:text-gray-700">← Back to Results</button>
          <div className="flex items-center gap-2 bg-teal-50 px-4 py-2 rounded-full">
            <span className="text-teal-500">✨</span>
            <span className="text-teal-600 font-medium text-sm">AI Assistant</span>
          </div>
        </div>

        {/* Messages */}
        <div className="flex-1 overflow-y-auto px-8 py-6 flex flex-col gap-4">
          {messages.map((msg, i) => (
            <div key={i} className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}>
              <div className={`max-w-lg px-5 py-4 rounded-2xl text-sm leading-relaxed ${
                msg.role === 'user'
                  ? 'bg-blue-900 text-white'
                  : 'bg-gray-100 text-gray-700'
              }`}>
                {msg.text}
              </div>
            </div>
          ))}
          {loading && (
            <div className="flex justify-start">
              <div className="bg-gray-100 px-5 py-4 rounded-2xl text-gray-400 text-sm animate-pulse">
                Analyzing...
              </div>
            </div>
          )}
          <div ref={bottomRef} />
        </div>

        {/* Suggestions */}
        <div className="px-8 pb-3 flex flex-wrap gap-2">
          {suggestions.map((s, i) => (
            <button key={i} onClick={() => sendMessage(s)}
              className="border border-gray-200 text-gray-600 text-xs px-4 py-2 rounded-full hover:bg-gray-50">
              {s}
            </button>
          ))}
        </div>

        {/* Input */}
        <div className="px-8 pb-6 flex items-center gap-3">
          <div className="flex-1 flex items-center border border-gray-200 rounded-2xl px-5 py-3 gap-3">
            <input
              className="flex-1 outline-none text-gray-700 text-sm"
              placeholder="Ask me anything about your DNA analysis..."
              value={input}
              onChange={e => setInput(e.target.value)}
              onKeyDown={e => e.key === 'Enter' && sendMessage()}
            />
            <button onClick={startListening}
              className={`text-lg transition-colors ${listening ? 'text-red-500 animate-pulse' : 'text-gray-400 hover:text-teal-500'}`}>
              🎤
            </button>
          </div>
          <button onClick={() => sendMessage()}
            className="bg-blue-900 text-white w-12 h-12 rounded-2xl flex items-center justify-center hover:bg-blue-800">
            ➤
          </button>
        </div>
      </div>

      {/* Sidebar */}
      <div className="w-72 border-l border-gray-100 p-6 flex flex-col gap-4">
        <h3 className="font-bold text-gray-900">Analysis Summary</h3>
        {[
          { icon: "⚡", label: "Sequence Type", value: "Human DNA" },
          { icon: "✓", label: "Mutation Status", value: "No mutations" },
          { icon: "↗", label: "Confidence", value: "98.7%" },
        ].map((item, i) => (
          <div key={i} className="flex items-center gap-3 border border-gray-100 rounded-xl p-3">
            <div className="w-9 h-9 bg-teal-50 rounded-lg flex items-center justify-center">{item.icon}</div>
            <div>
              <p className="text-gray-400 text-xs">{item.label}</p>
              <p className="font-bold text-gray-900 text-sm">{item.value}</p>
            </div>
          </div>
        ))}

        <div className="bg-gradient-to-br from-blue-900 to-teal-500 rounded-2xl p-4 mt-2">
          <h4 className="font-bold text-white mb-3">Quick Tips</h4>
          {[
            "Ask specific questions for detailed answers",
            "Request explanations in simpler terms",
            "Inquire about next steps or recommendations",
          ].map((tip, i) => (
            <p key={i} className="text-white/70 text-xs mb-2">• {tip}</p>
          ))}
        </div>
      </div>
    </div>
  )
}