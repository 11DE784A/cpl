#include "cpl_includes.h"
#include "cpl_defines.h"

#include "cpl_commons.h"
#include "cpl_format.h"

static const char *TUNITS[4] = {"s", "ms", "μs", "ns"};
static const char *BUNITS[5] = {"B", "KiB", "MiB", "GiB", "TiB"};

void strfsecs(char *str, double secs) {
	cpl_check(secs >= 0, "Negative times not supported");

	if (secs < 1) {
		int power = abs(log10(secs)) / 3;
		float fsecs = secs * pow(1000, power);
		sprintf(str, "%.4g %s", fsecs, TUNITS[power]);
	} else {
		sprintf(str, "%.4g %s", secs, TUNITS[0]);
	}
}

void strfbytes(char *str, int bytes) {
	cpl_check(bytes >= 0, "Number of bytes should be positive");

	int power = bytes > 0 ? log2(bytes) / 10 : 0;
	float fbytes = bytes / pow(1024, power);
	sprintf(str, "%.4g %s", fbytes, BUNITS[power]);
}

