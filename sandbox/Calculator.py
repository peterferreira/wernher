# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run -i '../Common.ipynb'

# <codecell>

μ = 398600
x0 = np.array([7000,-12124])
v0 = np.array([2.6679,4.6210])
r0 = np.sqrt(sum(r0**2))
s0 = np.sqrt(sum(v0**2))
vr0 = sum(x0 * v0) / r0
a = 1 / ((2 / r0) - ((s0**2) / μ))
α = 1 / a

vt0 = np.sqrt(s0**2 - vr0**2)
h = vt0 * r0
e = np.sqrt(1 + (r0 * vt0**2 / μ) * ((r0 * s0**2 / μ) - 2))

aa = h**2 / (μ * (1 - e**2))
αα = 1 / aa

print('''\
r0 = {r0:.0f}
s0 = {s0:.4f}
vr0 = {vr0:.4f}
a = {a:.0f}
α = {α:.4g}

vt0 = {vt0:.4f}
h = {h:.0f}
e = {e:.4f}
aa = {aa:.0f}
αα = {αα:.4g}
'''.format(
    r0 = r0,s0 = s0,vr0 = vr0,a = a,α = α,vt0 = vt0,h = h,e = e,
    aa = aa,αα = αα
))

