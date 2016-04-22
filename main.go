// Package probability provides a probability-theory toolbox.
package probability

// Sampler is a probability distribution capable of sampling.
type Sampler interface {
	Sample(Generator) float64
}

// Cumulator is a probability distribution capable of evaluating its CDF.
type Cumulator interface {
	Cumulate(float64) float64
}

// Decumulator is a probability distribution capable of evaluating the inverse
// of its CDF.
type Decumulator interface {
	Decumulate(float64) float64
}

// Distribution is a probability distribution.
type Distribution interface {
	Cumulator
	Decumulator
	Sampler
}

// Sample draws samples from the given distribution.
func Sample(distribution Sampler, generator Generator, count uint) []float64 {
	result := make([]float64, count)
	for i := range result {
		result[i] = distribution.Sample(generator)
	}
	return result
}

// Cumulate evaluates a CDF at a set of points.
func Cumulate(distribution Cumulator, x []float64) []float64 {
	result := make([]float64, len(x))
	for i := range result {
		result[i] = distribution.Cumulate(x[i])
	}
	return result
}

// Decumulate evaluates the inverse of a CDF at a set of points.
func Decumulate(distribution Decumulator, x []float64) []float64 {
	result := make([]float64, len(x))
	for i := range result {
		result[i] = distribution.Decumulate(x[i])
	}
	return result
}
