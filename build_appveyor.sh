
set -e
set -u

GENERATOR="$1"
AOPTION="$2"
BUILDDIR="$3"

(\
mkdir -p ${BUILDDIR} && \
cd ${BUILDDIR} && \
cmake -G "${GENERATOR}" -A "${AOPTION}" ../umfpack && \
cmake --build . --config Release -- //m //nologo //verbosity:minimal \
)

