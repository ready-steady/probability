// Package prob provides functions for working with probability distributions.
package prob

// Sampler represents a probability distribution capable of sampling.
type Sampler interface {
	Sample() float64
}

// Cumulator represents a probability distribution capable of evaluating its
// cumulative distribution function.
type Cumulator interface {
	CDF(float64) float64
}

// Inverter represents a probability distribution capable of evaluating the
// inverse of its cumulative distribution function.
type Inverter interface {
	InvCDF(float64) float64
}

// Distribution represents a probability distribution.
type Distribution interface {
	Cumulator
	Inverter
	Sampler
}

// Sample draws samples from the given distribution.
func Sample(dist Sampler, count uint32) []float64 {
	result := make([]float64, count)

	for i := range result {
		result[i] = dist.Sample()
	}

	return result
}

// CDF evaluates the cumulative distribution function of the given
// distribution.
func CDF(dist Cumulator, x []float64) []float64 {
	result := make([]float64, len(x))

	for i := range result {
		result[i] = dist.CDF(x[i])
	}

	return result
}

// InvCDF evaluates the inverse of the cumulative distribution function of the
// given distribution.
func InvCDF(dist Inverter, x []float64) []float64 {
	result := make([]float64, len(x))

	for i := range result {
		result[i] = dist.InvCDF(x[i])
	}

	return result
}
