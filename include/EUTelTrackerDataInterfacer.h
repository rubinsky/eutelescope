// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELTRACKERDATAINTERFACER_H
#define EUTELTRACKERDATAINTERFACER_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelGeometricPixel.h"


#ifdef USE_MARLIN
// marling includes ".h"
#include <marlin/Exceptions.h>
#endif

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <iomanip>

// template implementation
#include "EUTelTrackerDataInterfacer.hcc"
#include "EUTelTrackerDataInterfacer.tcc"

#endif

