// Package beta provides algorithms for working with beta distributions.
//
// https://en.wikipedia.org/wiki/Beta_distribution
package beta

import (
	"math"
)

// Self represents a particular distribution from the family.
type Self struct {
	α float64
	β float64
	a float64
	b float64
}

// New returns a beta distribution with α and β on [a, b].
func New(α, β, a, b float64) *Self {
	return &Self{α, β, a, b}
}

// CDF evaluates the CDF of the distribution.
func (s *Self) CDF(points []float64) []float64 {
	values := make([]float64, len(points))

	α, β, k, b := s.α, s.β, s.b-s.a, s.a
	logBeta := logBeta(α, β)

	for i, x := range points {
		values[i] = incBeta((x-b)/k, α, β, logBeta)
	}

	return values
}

// InvCDF evaluates the inverse CDF of the distribution.
func (s *Self) InvCDF(points []float64) []float64 {
	values := make([]float64, len(points))

	α, β, k, b := s.α, s.β, s.b-s.a, s.a
	logBeta := logBeta(α, β)

	for i, x := range points {
		values[i] = k*invIncBeta(x, α, β, logBeta) + b
	}

	return values
}

func incBeta(x, α, β, logBeta float64) float64 {
	// Author: John Burkardt
	// Source: http://people.sc.fsu.edu/~jburkardt/c_src/asa063/asa063.html
	// Modified: October 31, 2010

	const (
		acu = 0.1e-14
	)

	if x <= 0 {
		return 0
	}
	if 1 <= x {
		return 1
	}

	sum := α + β
	αx, βx := x, 1-x

	// Change the tail if necessary.
	var flip bool
	if α < sum*x {
		α, αx, β, βx = β, βx, α, αx
		flip = true
	}

	// Use Soper’s reduction formula.
	rx := αx / βx

	ns := int(β + βx*sum)
	if ns == 0 {
		rx = αx
	}

	ai := 1
	temp := β - float64(ai)
	term := 1.0

	value := 1.0

	for {
		term = term * temp * rx / (α + float64(ai))

		value += term

		temp = math.Abs(term)
		if temp <= acu && temp <= acu*value {
			break
		}

		ai++
		ns--

		if 0 < ns {
			temp = β - float64(ai)
		} else if ns == 0 {
			temp = β - float64(ai)
			rx = αx
		} else {
			temp = sum
			sum += 1
		}
	}

	value = value * math.Exp(α*math.Log(αx)+(β-1)*math.Log(βx)-logBeta) / α

	if flip {
		return 1 - value
	} else {
		return value
	}
}

func logBeta(x, y float64) float64 {
	z, _ := math.Lgamma(x + y)
	x, _ = math.Lgamma(x)
	y, _ = math.Lgamma(y)

	return x + y - z
}

func invIncBeta(x, α, β, logBeta float64) float64 {
	// Author: John Burkardt
	// Source: http://people.sc.fsu.edu/~jburkardt/c_src/asa109/asa109.html
	// Modified: September 17, 2014

	const (
		sae = -30
	)

	if x == 0 {
		return 0
	}
	if x == 1 {
		return 1
	}

	var iex int
	var acu, adj, fpu, g, h float64
	var prev, r, s, sq, t, tx, w, xin, y, yprev float64

	fpu = math.Pow10(sae)

	var flip bool
	if 0.5 < x {
		x = 1 - x
		α, β = β, α
		flip = true
	}

	var value float64

	// Calculate the initial approximation.
	r = math.Sqrt(-math.Log(x * x))
	y = r - (2.30753+0.27061*r)/(1+(0.99229+0.04481*r)*r)

	if 1 < α && 1 < β {
		r = (y*y - 3) / 6
		s = 1 / (2*α - 1)
		t = 1 / (2*β - 1)
		h = 2 / (s + t)
		w = y*math.Sqrt(h+r)/h - (t-s)*(r+5/6-2/(3*h))
		value = α / (α + β*math.Exp(w+w))
	} else {
		r = β + β
		t = 1 / (9 * β)
		t = r * math.Pow(1-t+y*math.Sqrt(t), 3)

		if t <= 0 {
			value = 1 - math.Exp((math.Log((1-x)*β)+logBeta)/β)
		} else {
			t = (4*α + r - 2) / t

			if t <= 1 {
				value = math.Exp((math.Log(x*α) + logBeta) / α)
			} else {
				value = 1 - 2/(t+1)
			}
		}
	}

	// Solve by a modified Newton–Raphson method.
	r = 1 - α
	t = 1 - β
	yprev = 0
	sq = 1
	prev = 1

	if value < 0.0001 {
		value = 0.0001
	}

	if 0.9999 < value {
		value = 0.9999
	}

	iex = int(-5/α/α - 1/math.Pow(x, 0.2) - 13)
	if iex < sae {
		iex = sae
	}

	acu = math.Pow10(iex)

outer:
	for {
		y = incBeta(value, α, β, logBeta)

		xin = value
		y = (y - x) * math.Exp(logBeta+r*math.Log(xin)+t*math.Log(1-xin))

		if y*yprev <= 0 {
			if sq > fpu {
				prev = sq
			} else {
				prev = fpu
			}
		}

		g = 1

		for {
			for {
				adj = g * y
				sq = adj * adj

				if sq < prev {
					tx = value - adj

					if 0 <= tx && tx <= 1 {
						break
					}
				}
				g = g / 3
			}

			if prev <= acu || y*y <= acu {
				break outer
			}

			if tx != 0 && tx != 1 {
				break
			}

			g = g / 3
		}

		if tx == value {
			break
		}

		value = tx
		yprev = y
	}

	if flip {
		return 1 - value
	} else {
		return value
	}
}
