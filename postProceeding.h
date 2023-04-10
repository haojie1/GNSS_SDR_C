#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector_int.h>
#include "track.h"
#include "config.h"
#include "initSetting.h"
#include "findPreamble.h"
#include "ephemeris.h"
#include "satpos.h"
#include "common.h"
#include "cart2geo.h"
struct posLLH getPos(struct settings receiverSetting) ;
