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

	for i, x := range points {
		values[i] = ibeta(α, β, (x-b)/k)
	}

	return values
}

// InvCDF evaluates the inverse CDF of the distribution.
func (s *Self) InvCDF(points []float64) []float64 {
	values := make([]float64, len(points))

	α, β, k, b := s.α, s.β, s.b-s.a, s.a

	for i, p := range points {
		values[i] = k*invbeta(p, α, β) + b
	}

	return values
}

func ibeta(α, β, x float64) float64 {
	// Author: John Burkardt
	// Source: http://people.sc.fsu.edu/~jburkardt/c_src/prob/prob.html

	const (
		max = 1000
		tol = 1e-07
	)

	if x <= 0 {
		return 0
	} else if 1 <= x {
		return 1
	}

	αx, βx := x, 1-x

	var flip bool
	if α < (α+β)*x {
		α, αx, β, βx = β, βx, α, αx
		flip = true
	}

	// https://en.wikipedia.org/wiki/Integration_by_reduction_formulae
	rx := αx
	n := int(β + βx*(α+β))
	if n != 0 {
		rx /= βx
	}

	value := 1.0

	for i, temp, psq, term := 1, β-1, α+β, 1.0; ; {
		term = term * temp * rx / (α + float64(i))

		value += term
		temp = math.Abs(term)

		if temp <= tol && temp <= tol*value {
			break
		}

		i++

		if i > max {
			// FIXME: Might not converge. Panic?
			break
		}

		n--

		if 0 <= n {
			temp = β - float64(i)
			if n == 0 {
				rx = αx
			}
		} else {
			temp = psq
			psq += 1
		}
	}

	// http://dlmf.nist.gov/8.17#v
	value = value * math.Exp(α*math.Log(αx)+(β-1)*math.Log(βx)) / (beta(α, β) * α)

	if flip {
		return 1 - value
	} else {
		return value
	}
}

func beta(x, y float64) float64 {
	z, _ := math.Lgamma(x + y)
	x, _ = math.Lgamma(x)
	y, _ = math.Lgamma(y)

	return math.Exp(x + y - z)
}

func invbeta(cdf, p, q float64) float64 {
	const (
		sae = -37
	)

	var flip bool
	var iex int
	var a, acu, adj, beta_log, fpu, g, h float64
	var pp, prev, qq, r, s, sq, t, tx, value, w, xin, y, yprev float64

	fpu = math.Pow10(sae)

	if cdf == 0 {
		return 0
	}
	if cdf == 1 {
		return 1
	}

	if 0.5 < cdf {
		a = 1 - cdf
		pp = q
		qq = p
		flip = true
	} else {
		a = cdf
		pp = p
		qq = q
	}

	// Calculate the initial approximation.
	r = math.Sqrt(-math.Log(a * a))

	y = r - (2.30753+0.27061*r)/(1+(0.99229+0.04481*r)*r)

	if 1 < pp && 1 < qq {
		r = (y*y - 3) / 6
		s = 1 / (pp + pp - 1)
		t = 1 / (qq + qq - 1)
		h = 2 / (s + t)
		w = y*math.Sqrt(h+r)/h - (t-s)*(r+5/6-2/(3*h))
		value = pp / (pp + qq*math.Exp(w+w))
	} else {
		r = qq + qq
		t = 1 / (9 * qq)
		t = r * math.Pow(1-t+y*math.Sqrt(t), 3)

		if t <= 0 {
			value = 1 - math.Exp((math.Log((1-a)*qq)+beta_log)/qq)
		} else {
			t = (4*pp + r - 2) / t

			if t <= 1 {
				value = math.Exp((math.Log(a*pp) + beta_log) / pp)
			} else {
				value = 1 - 2/(t+1)
			}
		}
	}

	// A modified Newton–Raphson method
	r = 1 - pp
	t = 1 - qq
	yprev = 0
	sq = 1
	prev = 1

	if value < 0.0001 {
		value = 0.0001
	}

	if 0.9999 < value {
		value = 0.9999
	}

	iex = int(-5/pp/pp - 1/math.Pow(a, 0.2) - 13)
	if sae > iex {
		iex = sae
	}

	acu = math.Pow10(iex)

	for {
		y = ibeta(pp, qq, value)

		xin = value
		y = (y - a) * math.Exp(beta_log+r*math.Log(xin)+t*math.Log(1-xin))

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

			if prev <= acu {
				if flip {
					value = 1 - value
				}
				return value
			}

			if y*y <= acu {
				if flip {
					value = 1 - value
				}
				return value
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
		value = 1 - value
	}

	return value
}
