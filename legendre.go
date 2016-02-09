package quadrature

import (
	"math"
)

// Legendre computes the Gauss–Legendre rule.
func Legendre(order uint) (x, w []float64) {
	x = make([]float64, order)
	w = make([]float64, order)

	e1 := float64(order * (order + 1))

	m := (order + 1) / 2

	for i := uint(1); i <= m; i++ {
		mp1mi := m + 1 - i

		t := float64(4*i-1) * math.Pi / float64(4*order+2)

		x0 := math.Cos(t) * (1.0 - (1.0-1.0/float64(order))/float64(8*order*order))

		pkm1 := 1.0
		pk := x0

		for k := uint(2); k <= order; k++ {
			pkp1 := 2.0*x0*pk - pkm1 - (x0*pk-pkm1)/float64(k)
			pkm1 = pk
			pk = pkp1
		}

		d1 := float64(order) * (pkm1 - x0*pk)

		dpn := d1 / (1.0 - x0*x0)
		d2pn := (2.0*x0*dpn - e1*pk) / (1.0 - x0*x0)
		d3pn := (4.0*x0*d2pn + (2.0-e1)*dpn) / (1.0 - x0*x0)
		d4pn := (6.0*x0*d3pn + (6.0-e1)*d2pn) / (1.0 - x0*x0)

		u := pk / dpn
		v := d2pn / dpn

		h := -u * (1.0 + 0.5*u*(v+u*(v*v-d3pn/(3.0*dpn))))

		p := pk + h*(dpn+0.5*h*(d2pn+h/3.0*(d3pn+0.25*h*d4pn)))
		dp := dpn + h*(d2pn+0.5*h*(d3pn+h*d4pn/3.0))
		h = h - p/dp

		xtemp := x0 + h
		x[mp1mi-1] = xtemp

		fx := d1 - h*e1*(pk+0.5*h*(dpn+h/3.0*(d2pn+0.25*h*(d3pn+0.2*h*d4pn))))

		w[mp1mi-1] = 2.0 * (1.0 - xtemp*xtemp) / (fx * fx)
	}

	if order%2 == 1 {
		x[0] = 0.0
	}

	nmove := (order + 1) / 2
	ncopy := order - nmove
	for i := uint(1); i <= nmove; i++ {
		iback := order + 1 - i
		x[iback-1] = x[iback-ncopy-1]
		w[iback-1] = w[iback-ncopy-1]
	}

	for i := uint(1); i <= order-nmove; i++ {
		x[i-1] = -x[order-i]
		w[i-1] = w[order-i]
	}

	return
}
