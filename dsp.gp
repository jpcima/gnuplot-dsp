sinc(x) = sin(x) / x
normsinc(x) = sinc(pi * x)

sincm(x, m) = normsinc(x) / (m * normsinc(x / m))

hav(x) = .5 * (1. - cos(x))

fmod(x, a) = x - a * int(x / a)
logbase(x, b) = log(x) / log(b)

amp2db(x) = 20. * logbase(x, 10.)
db2amp(x) = 10. ** (x / 20.)

notefreq(x) = 440. * 2.**((x - 69.) / 12.)
freqnote(x) = 12. * logbase(x / 440., 2.) + 69.

besI0(x) = (abs(x) < 3.75) ? (1. +                                                  \
                              (((x / 3.75)**2.)**1. * 3.5156229) +                  \
                              (((x / 3.75)**2.)**2. * 3.0899424) +                  \
                              (((x / 3.75)**2.)**3. * 1.2067492) +                  \
                              (((x / 3.75)**2.)**4. * 0.2659732) +                  \
                              (((x / 3.75)**2.)**5. * 0.360768e-1) +                \
                              (((x / 3.75)**2.)**6. * 0.45813e-2))                  \
           : ((exp(abs(x)) / sqrt(abs(x))) * (0.39894228 +                            \
                                              ((3.75 / abs(x))**1. * 0.1328592e-1) +  \
                                              ((3.75 / abs(x))**2. * 0.225319e-2) +   \
                                              ((3.75 / abs(x))**3. * -0.157565e-2) +  \
                                              ((3.75 / abs(x))**4. * 0.916281e-2) +   \
                                              ((3.75 / abs(x))**5. * -0.2057706e-1) + \
                                              ((3.75 / abs(x))**6. * 0.2635537e-1) +  \
                                              ((3.75 / abs(x))**7. * -0.1647633e-1) + \
                                              ((3.75 / abs(x))**8. * 0.392377e-2)))

window_hann(x) = (x < 0. || x > 1) ? 0.           \
               : (.5 * (1. - cos(2. * pi * x)))

window_tukey(x, a) = (x < 0. || x > 1) ? 0.                                                \
                   : (x < a / 2.) ? ((1. + cos(pi * (2. * x / a - 1.))) / 2.)              \
                   : (x > 1. - a / 2.) ? ((1. + cos(pi * (2. * (x - 1.) / a + 1.))) / 2.)  \
                   : 1

window_kaiser(x, a) = (x < 0. || x > 1) ? 0.                             \
                    : (besI0((a * pi) * sqrt(1. - (2. * x - 1.)**2.)) /  \
                       besI0(pi * a))

window_triangular(x) = (x < 0. || x > 1) ? 0.  \
                     : (x < .5) ? x            \
                     : (1. - x)

window_hann(x) = (x < 0. || x > 1) ? 0.  \
               : hav(2. * pi * x)

window_blackman_(x, a) = (x < 0. || x > 1) ? 0.  \
                       : (1. - a) / 2. - .5 * cos(2 * pi * x) + (a / 2.) * cos(4.0 * pi * x)
window_blackman(x) = window_blackman_(x, 0.16)

window_cosine(x) = (x < 0. || x > 1) ? 0.  \
                 : sin(pi * x)

window_lanczos(x) = (x < 0. || x > 1) ? 0.  \
                  : normsinc(2. * x - 1.)

periodicity_(x, f) = fmod(x, 1. / f) / (1. / f)
periodicity(x, f) = (periodicity_(x, f) < 0.) ? (periodicity_(x, f) + 1.0) : periodicity_(x, f)

wave_square(x, f) = (periodicity(x, f) < .5) ? -1. : 1.
wave_triangle(x, f) = 4. * ((periodicity(x, f) < .5) ? periodicity(x, f) : (1. - periodicity(x, f))) - 1.
wave_ramp(x, f) = 2. * periodicity(x, f) - 1.
wave_saw(x, f) = 2. * (1. - periodicity(x, f)) - 1.
wave_sine(x, f) = sin(2. * pi * f * x)
