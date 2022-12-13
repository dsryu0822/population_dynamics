const a₁ = 15
const a₂ = 20
const x = (-a₁):0.1:a₂
const x30 = x.+30

Pearsonian1(x, y₀, m₁, m₂) = y₀*((1 + (x/a₁))^m₁)*((1 - (x/a₂))^m₂)