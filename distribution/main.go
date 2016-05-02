// Package provides probability distributions.
package distribution

// Continuous is a continuous distribution.
type Continuous interface {
	Cumulator
	Denser
	Inverter
	Sampler
}

// Cumulator is a distribution capable of evaluating its CDF.
type Cumulator interface {
	Cumulate(float64) float64
}

// Denser is a distribution capable of evaluating its PDF.
type Denser interface {
	Dense(float64) float64
}

// Inverter is a distribution capable of evaluating the inverse of its CDF.
type Inverter interface {
	Invert(float64) float64
}

// Sampler is a distribution capable of sampling.
type Sampler interface {
	Sample(Generator) float64
}

// Cumulate evaluates a CDF at a set of points.
func Cumulate(distribution Cumulator, x []float64) []float64 {
	result := make([]float64, len(x))
	for i := range result {
		result[i] = distribution.Cumulate(x[i])
	}
	return result
}

// Dense evaluates a PDF at a set of points.
func Dense(distribution Denser, x []float64) []float64 {
	result := make([]float64, len(x))
	for i := range result {
		result[i] = distribution.Dense(x[i])
	}
	return result
}

// Invert evaluates the inverse of a CDF at a set of points.
func Invert(distribution Inverter, x []float64) []float64 {
	result := make([]float64, len(x))
	for i := range result {
		result[i] = distribution.Invert(x[i])
	}
	return result
}

// Sample draws samples from the given distribution.
func Sample(distribution Sampler, generator Generator, count uint) []float64 {
	result := make([]float64, count)
	for i := range result {
		result[i] = distribution.Sample(generator)
	}
	return result
}
