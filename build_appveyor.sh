
set -e
set -u

GENERATOR="$1"
AOPTION="$2"
TOOLSET="$3"
BUILDDIR="$4"

(\
mkdir -p ${BUILDDIR} && \
cd ${BUILDDIR} && \
cmake -G "${GENERATOR}" -A "${AOPTION}" -T "${TOOLSET}" ../umfpack && \
cmake --build . --config Release --parallel 4 -- //nologo //verbosity:minimal \
)

