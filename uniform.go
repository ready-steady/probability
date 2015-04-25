package probability

// Uniform represents a uniform distribution.
type Uniform struct {
	a float64
	b float64
}

// NewUniform returns a uniform distribution on [a, b].
func NewUniform(a, b float64) *Uniform {
	return &Uniform{a, b}
}

// Sample draws a sample from the distribution.
func (u *Uniform) Sample(generator Generator) float64 {
	return (u.b-u.a)*generator.Float64() + u.a
}
