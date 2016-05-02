package distribution

import (
	"math"
)

var (
	sqrt2Pi = math.Sqrt(2.0 * math.Pi)
)

// Gaussian is a Gaussian distribution.
type Gaussian struct {
	μ float64
	σ float64
}

// NewGaussian returns a Gaussian distribution with mean μ and standard
// deviation σ.
func NewGaussian(μ, σ float64) *Gaussian {
	return &Gaussian{μ, σ}
}

// Cumulate evaluates the CDF.
func (self *Gaussian) Cumulate(x float64) float64 {
	return (1.0 + math.Erf((x-self.μ)/(self.σ*math.Sqrt2))) / 2.0
}

// Dense evaluates the PDF.
func (self *Gaussian) Dense(x float64) float64 {
	μ, σ := self.μ, self.σ
	return math.Exp(-(x-μ)*(x-μ)/(2.0*σ*σ)) / (sqrt2Pi * σ)
}

// Invert evaluates the inverse of the CDF.
//
// The code is based on a C implementation by John Burkardt.
// http://people.sc.fsu.edu/~jburkardt/c_src/asa241/asa241.html
func (self *Gaussian) Invert(p float64) (x float64) {
	const (
		const1 = 0.180625
		const2 = 1.6
		split1 = 0.425
		split2 = 5.0
	)

	if p <= 0.0 {
		return math.Inf(-1.0)
	}
	if 1.0 <= p {
		return math.Inf(1.0)
	}

	q := p - 0.5

	if math.Abs(q) <= split1 {
		x = const1 - q*q
		x = self.μ + self.σ*q*poly7(coefA, x)/poly7(coeffB, x)
		return
	}

	if q < 0.0 {
		x = p
	} else {
		x = 1.0 - p
	}

	x = math.Sqrt(-math.Log(x))

	if x <= split2 {
		x -= const2
		x = poly7(coefC, x) / poly7(coefD, x)
	} else {
		x -= split2
		x = poly7(coefE, x) / poly7(coefF, x)
	}

	if q < 0.0 {
		x = self.μ - self.σ*x
	} else {
		x = self.μ + self.σ*x
	}

	return
}

// Sample draws a sample.
func (self *Gaussian) Sample(generator Generator) float64 {
	return self.μ + self.σ*generator.NormFloat64()
}

func poly7(coef []float64, x float64) (y float64) {
	for i := 8 - 1; 0 <= i; i-- {
		y = y*x + coef[i]
	}
	return
}

var coefA = []float64{
	3.3871328727963666080,
	1.3314166789178437745e+2,
	1.9715909503065514427e+3,
	1.3731693765509461125e+4,
	4.5921953931549871457e+4,
	6.7265770927008700853e+4,
	3.3430575583588128105e+4,
	2.5090809287301226727e+3,
}

var coeffB = []float64{
	1.0,
	4.2313330701600911252e+1,
	6.8718700749205790830e+2,
	5.3941960214247511077e+3,
	2.1213794301586595867e+4,
	3.9307895800092710610e+4,
	2.8729085735721942674e+4,
	5.2264952788528545610e+3,
}

var coefC = []float64{
	1.42343711074968357734,
	4.63033784615654529590,
	5.76949722146069140550,
	3.64784832476320460504,
	1.27045825245236838258,
	2.41780725177450611770e-1,
	2.27238449892691845833e-2,
	7.74545014278341407640e-4,
}

var coefD = []float64{
	1.0,
	2.05319162663775882187,
	1.67638483018380384940,
	6.89767334985100004550e-1,
	1.48103976427480074590e-1,
	1.51986665636164571966e-2,
	5.47593808499534494600e-4,
	1.05075007164441684324e-9,
}

var coefE = []float64{
	6.65790464350110377720,
	5.46378491116411436990,
	1.78482653991729133580,
	2.96560571828504891230e-1,
	2.65321895265761230930e-2,
	1.24266094738807843860e-3,
	2.71155556874348757815e-5,
	2.01033439929228813265e-7,
}

var coefF = []float64{
	1.0,
	5.99832206555887937690e-1,
	1.36929880922735805310e-1,
	1.48753612908506148525e-2,
	7.86869131145613259100e-4,
	1.84631831751005468180e-5,
	1.42151175831644588870e-7,
	2.04426310338993978564e-15,
}
