#pragma once

#include <astray/benchmark/benchmark.hpp>

#include <astray/core/geodesic.hpp>
#include <astray/core/metric.hpp>
#include <astray/core/observer.hpp>
#include <astray/core/ray_tracer.hpp>

#include <astray/math/ode/ode.hpp>
#include <astray/math/angle.hpp>
#include <astray/math/constants.hpp>
#include <astray/math/coordinate_system.hpp>
#include <astray/math/linear_algebra.hpp>

#include <astray/media/image.hpp>
#include <astray/media/video.hpp>

#include <astray/metrics/boyer_lindquist/kerr.hpp>
#include <astray/metrics/cartesian/alcubierre.hpp>
#include <astray/metrics/cartesian/bessel.hpp>
#include <astray/metrics/cartesian/de_sitter.hpp>
#include <astray/metrics/cartesian/goedel.hpp>
#include <astray/metrics/cartesian/kastor_traschen.hpp>
#include <astray/metrics/cartesian/minkowski.hpp>
#include <astray/metrics/cylindrical/einstein_rosen_weber_wheeler_bonnor.hpp>
#include <astray/metrics/spherical/barriola_vilenkin.hpp>
#include <astray/metrics/spherical/bertotti_kasner.hpp>
#include <astray/metrics/spherical/friedman_lemaitre_robertson_walker.hpp>
#include <astray/metrics/spherical/janis_newman_winicour.hpp>
#include <astray/metrics/spherical/kottler.hpp>
#include <astray/metrics/spherical/morris_thorne.hpp>
#include <astray/metrics/spherical/reissner_nordstroem.hpp>
#include <astray/metrics/spherical/reissner_nordstroem_extreme_dihole.hpp>
#include <astray/metrics/spherical/schwarzschild.hpp>
#include <astray/metrics/spherical/schwarzschild_cosmic_string.hpp>

#include <astray/parallel/distributed_device.hpp>
#include <astray/parallel/shared_device.hpp>