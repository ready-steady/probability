package beta

import (
	"testing"

	"github.com/go-math/support/assert"
)

func TestCDF(t *testing.T) {
	points := []float64{
		-1.00, -0.85, -0.70, -0.55, -0.40, -0.25, -0.10, 0.05, 0.20, 0.35,
		0.50, 0.65, 0.80, 0.95, 1.10, 1.25, 1.40, 1.55, 1.70, 1.85, 2.00,
	}

	values := []float64{
		0.000000000000000e+00,
		1.401875000000000e-02,
		5.230000000000002e-02,
		1.095187500000000e-01,
		1.807999999999999e-01,
		2.617187500000001e-01,
		3.483000000000000e-01,
		4.370187500000001e-01,
		5.248000000000003e-01,
		6.090187500000001e-01,
		6.875000000000000e-01,
		7.585187500000001e-01,
		8.208000000000000e-01,
		8.735187499999999e-01,
		9.163000000000000e-01,
		9.492187500000000e-01,
		9.728000000000000e-01,
		9.880187500000001e-01,
		9.963000000000000e-01,
		9.995187500000000e-01,
		1.000000000000000e+00,
	}

	assert.AlmostEqual(New(2, 3, -1, 2).CDF(points), values, t)
}

func TestInvCDF(t *testing.T) {
	points := []float64{
		0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
		0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
	}

	values := []float64{
		3.000000000000000e+00,
		3.025320565519104e+00,
		3.051316701949486e+00,
		3.078045554270711e+00,
		3.105572809000084e+00,
		3.133974596215561e+00,
		3.163339973465924e+00,
		3.193774225170145e+00,
		3.225403330758517e+00,
		3.258380151290432e+00,
		3.292893218813452e+00,
		3.329179606750063e+00,
		3.367544467966324e+00,
		3.408392021690038e+00,
		3.452277442494834e+00,
		3.500000000000000e+00,
		3.552786404500042e+00,
		3.612701665379257e+00,
		3.683772233983162e+00,
		3.776393202250021e+00,
		4.000000000000000e+00,
	}

	assert.AlmostEqual(New(1, 2, 3, 4).InvCDF(points), values, t)
}
