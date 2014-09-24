// Package prob provides functions for working with probability distributions.
package prob

// Sampler represents a probability distribution capable of sampling.
type Sampler interface {
	Sample() float64
}

// CDFer represents a probability distribution capable of evaluating its CDF.
type CDFer interface {
	CDF(float64) float64
}

// InvCDFer represents a probability distribution capable of evaluating the
// inverse of its CDF.
type InvCDFer interface {
	InvCDF(float64) float64
}

// Sample draws samples from the given distribution.
func Sample(dist Sampler, count uint32) []float64 {
	result := make([]float64, count)

	for i := range result {
		result[i] = dist.Sample()
	}

	return result
}

// CDF evalues the CDF of the given distribution.
func CDF(dist CDFer, x[]float64) []float64 {
	result := make([]float64, len(x))

	for i := range result {
		result[i] = dist.CDF(x[i])
	}

	return result
}

// InvCDF evalues the inverse of the CDF of the given distribution.
func InvCDF(dist InvCDFer, x[]float64) []float64 {
	result := make([]float64, len(x))

	for i := range result {
		result[i] = dist.InvCDF(x[i])
	}

	return result
}
