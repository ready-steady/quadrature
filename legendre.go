package quadrature

import (
	"math"
)

// Legendre computes the Gauss–Legendre rule.
func Legendre(order uint, a, b float64) (x, w []float64) {
	x, w = legendre(order)
	rescale(order, x, w, a, b)
	return
}

func legendre(order uint) (x, w []float64) {
	x = make([]float64, order)
	w = make([]float64, order)

	e := float64(order * (order + 1))
	m := (order + 1) / 2
	for i := uint(1); i <= m; i++ {
		x0 := math.Cos(float64(4*i-1)*math.Pi/float64(4*order+2)) *
			(1.0 - (1.0-1.0/float64(order))/float64(8*order*order))

		p, q := x0, 1.0
		for j := uint(2); j <= order; j++ {
			p, q = 2.0*x0*p-q-(x0*p-q)/float64(j), p
		}

		d0 := float64(order) * (q - x0*p)
		d1 := d0 / (1.0 - x0*x0)
		d2 := (2.0*x0*d1 - e*p) / (1.0 - x0*x0)
		d3 := (4.0*x0*d2 + (2.0-e)*d1) / (1.0 - x0*x0)
		d4 := (6.0*x0*d3 + (6.0-e)*d2) / (1.0 - x0*x0)

		u, v := p/d1, d2/d1
		h := -u * (1.0 + 0.5*u*(v+u*(v*v-d3/(3.0*d1))))
		h -= (p + h*(d1+0.5*h*(d2+h/3.0*(d3+0.25*h*d4)))) / (d1 + h*(d2+0.5*h*(d3+h*d4/3.0)))

		f := d0 - h*e*(p+0.5*h*(d1+h/3.0*(d2+0.25*h*(d3+0.2*h*d4))))

		x[m-i] = x0 + h
		w[m-i] = 2.0 * (1.0 - x[m-i]*x[m-i]) / (f * f)
	}

	if order%2 == 1 {
		x[0] = 0.0
	}

	for i := uint(1); i <= m; i++ {
		x[order-i] = x[m-i]
		w[order-i] = w[m-i]
	}
	for i := uint(1); i <= order-m; i++ {
		x[i-1] = -x[order-i]
		w[i-1] = w[order-i]
	}

	return
}

func rescale(order uint, x, w []float64, a, b float64) {
	for i := uint(0); i < order; i++ {
		x[i] = ((a + b) + (b-a)*x[i]) / 2.0
		w[i] = (b - a) * w[i] / 2.0
	}
}
