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
	logB := logBeta(α, β)

	for i, x := range points {
		values[i] = incBeta((x-b)/k, α, β, logB)
	}

	return values
}

// InvCDF evaluates the inverse CDF of the distribution.
func (s *Self) InvCDF(points []float64) []float64 {
	values := make([]float64, len(points))

	α, β, k, b := s.α, s.β, s.b-s.a, s.a
	logB := logBeta(α, β)

	for i, p := range points {
		values[i] = k*invIncBeta(p, α, β, logB) + b
	}

	return values
}

func incBeta(x, p, q, logB float64) float64 {
	// The code is based on a C implementation by John Burkardt.
	// http://people.sc.fsu.edu/~jburkardt/c_src/asa109/asa109.html

	const (
		acu = 0.1e-14
	)

	if x <= 0 {
		return 0
	}
	if 1 <= x {
		return 1
	}

	sum := p + q
	px, qx := x, 1-x

	// Change the tail if necessary.
	var flip bool
	if p < sum*x {
		p, px, q, qx = q, qx, p, px
		flip = true
	}

	// Use Soper’s reduction formula.
	rx := px / qx

	ns := int(q + qx*sum)
	if ns == 0 {
		rx = px
	}

	ai := 1
	temp := q - float64(ai)
	term := 1.0

	α := 1.0

	for {
		term = term * temp * rx / (p + float64(ai))

		α += term

		temp = math.Abs(term)
		if temp <= acu && temp <= acu*α {
			break
		}

		ai++
		ns--

		if 0 < ns {
			temp = q - float64(ai)
		} else if ns == 0 {
			temp = q - float64(ai)
			rx = px
		} else {
			temp = sum
			sum += 1
		}
	}

	// Applied Statistics. Algorithm AS 109
	// http://www.jstor.org/discover/10.2307/2346887
	α = α * math.Exp(p*math.Log(px)+(q-1)*math.Log(qx)-logB) / p

	if flip {
		return 1 - α
	} else {
		return α
	}
}

func invIncBeta(α, p, q, logB float64) float64 {
	// The code is based on a C implementation by John Burkardt.
	// http://people.sc.fsu.edu/~jburkardt/c_src/asa109/asa109.html

	const (
		// Applied Statistics. Algorithm AS R83
		// http://www.jstor.org/discover/10.2307/2347779
		//
		// The machine-dependent smallest allowable exponent of 10 to avoid
		// floating-point underflow error.
		sae = -30
	)

	if α <= 0 {
		return 0
	}
	if 1 <= α {
		return 1
	}

	var flip bool
	if 0.5 < α {
		α = 1 - α
		p, q = q, p
		flip = true
	}

	// An approximation x₀ to x if found from (cf. Scheffé and Tukey, 1944)
	//
	//     (1 + x₀)/(1 - x₀) = (4*p + 2*q - 2)/χ²(α)
	//
	// where χ²(α) is the upper α point of the χ² distribution with 2*q degrees
	// of freedom and is obtained from Wilson and Hilferty’s approximation (cf.
	// Wilson and Hilferty, 1931)
	//
	//     χ²(α) = 2*q*(1 - 1/(9*q) + y(α) * sqrt(1/(9*q)))**3,
	//
	// y(α) being Hastings’ approximation (cf. Hastings, 1955) for the upper α
	// point of the standard normal distribution. If χ²(α) < 0, then
	//
	//     x₀ = 1 - ((1 - α)*q*B(p, q))**(1/q).
	//
	// Again if (4*p + 2*q - 2)/χ²(α) does not exceed 1, x₀ is obtained from
	//
	//     x₀ = (α*p*B(p, q))**(1/p).
	//
	// Applied Statistics. Algorithm AS 46
	// http://www.jstor.org/discover/10.2307/2346798
	var x float64

	var y, r, t float64

	r = math.Sqrt(-math.Log(α * α))
	y = r - (2.30753+0.27061*r)/(1+(0.99229+0.04481*r)*r)

	if 1 < p && 1 < q {
		// For p and q > 1, the approximation given by Carter (1947), which
		// improves the Fisher–Cochran formula, is generally better. For other
		// values of p and q en empirical investigation has shown that the
		// approximation given in AS 64 is adequate.
		//
		// Applied Statistics. Algorithm AS 109
		// http://www.jstor.org/discover/10.2307/2346887
		r = (y*y - 3) / 6
		s := 1 / (2*p - 1)
		t = 1 / (2*q - 1)
		h := 2 / (s + t)
		w := y*math.Sqrt(h+r)/h - (t-s)*(r+5/6-2/(3*h))
		x = p / (p + q*math.Exp(2*w))
	} else {
		t = 1 / (9 * q)
		t = 2 * q * math.Pow(1-t+y*math.Sqrt(t), 3)
		if t <= 0 {
			x = 1 - math.Exp((math.Log((1-α)*q)+logB)/q)
		} else {
			t = 2 * (2*p + q - 1) / t
			if t <= 1 {
				x = math.Exp((math.Log(α*p) + logB) / p)
			} else {
				x = 1 - 2/(t+1)
			}
		}
	}

	if x < 0.0001 {
		x = 0.0001
	} else if 0.9999 < x {
		x = 0.9999
	}

	// The final solution is obtained by the Newton–Raphson method from the
	// relation
	//
	//     x[i] = x[i-1] - f(x[i-1])/f'(x[i-1])
	//
	// where
	//
	//     f(x) = I(x, p, q) - α.
	//
	// Applied Statistics. Algorithm AS 46
	// http://www.jstor.org/discover/10.2307/2346798
	r = 1 - p
	t = 1 - q
	yprev := 0.0
	sq := 1.0
	prev := 1.0

	// Applied Statistics. Algorithm AS R83
	// http://www.jstor.org/discover/10.2307/2347779
	fpu := math.Pow10(sae)
	acu := fpu
	if e := int(-5/p/p - 1/math.Pow(α, 0.2) - 13); e > sae {
		acu = math.Pow10(e)
	}

	var tx, g, adj float64

outer:
	for {
		// Applied Statistics. Algorithm AS 109
		// http://www.jstor.org/discover/10.2307/2346887
		y = incBeta(x, p, q, logB)
		y = (y - α) * math.Exp(logB+r*math.Log(x)+t*math.Log(1-x))

		// Applied Statistics. Algorithm AS R83
		// http://www.jstor.org/discover/10.2307/2347779
		if y*yprev <= 0 {
			prev = math.Max(sq, fpu)
		}

		g = 1

		for {
			for {
				adj = g * y
				sq = adj * adj

				if sq < prev {
					tx = x - adj

					if 0 <= tx && tx <= 1 {
						break
					}
				}
				g /= 3
			}

			if prev <= acu || y*y <= acu {
				x = tx
				break outer
			}

			if tx != 0 && tx != 1 {
				break
			}

			g /= 3
		}

		if tx == x {
			break
		}

		x = tx
		yprev = y
	}

	if flip {
		return 1 - x
	} else {
		return x
	}
}

func logBeta(x, y float64) float64 {
	z, _ := math.Lgamma(x + y)
	x, _ = math.Lgamma(x)
	y, _ = math.Lgamma(y)

	return x + y - z
}
