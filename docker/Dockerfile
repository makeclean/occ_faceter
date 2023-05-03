FROM ubuntu:latest AS ubuntu_with_occ

# install runtime packages (OCCT pulls in a lot of dependencies)
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates \
      libhdf5-103-1 \
      liblapack3 \
      libocct-data-exchange-7.5 \
      libocct-draw-7.5 \
      libocct-foundation-7.5 \
      libocct-modeling-data-7.5 \
      libtbb12 \
 && rm -rf /var/lib/apt/lists/*

# build code in larger inner image containing compilers and -dev packages
FROM ubuntu_with_occ AS inner

ENV DEBIAN_FRONTEND=noninteractive

# get system packages installed
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      automake autoconf libtool \
      build-essential \
      cmake \
      gfortran \
      git ca-certificates \
      libcgal-dev \
      libeigen3-dev \
      libhdf5-dev \
      liblapack-dev \
      libocct-data-exchange-dev \
      libocct-draw-dev \
      libocct-foundation-dev \
      libocct-modeling-data-dev \
      libtbb-dev \
      occt-misc \
      python3 \
      tzdata \
 && rm -rf /var/lib/apt/lists/*

# clone repos once to speed rebuilds up
RUN git clone --depth=1 https://bitbucket.org/fathomteam/moab.git /moab \
 && cd /moab \
 && autoreconf -fi

RUN git clone --depth=1 --recurse-submodules --shallow-submodules https://github.com/svalinn/DAGMC /DAGMC

RUN git clone --depth=1 --recurse-submodules --shallow-submodules https://github.com/ukaea/overlap_checker /overlap_checker

RUN git clone --depth=1 --branch=new-materials https://github.com/johnnonweiler/occ_faceter /occ_faceter

# build and install MOAB
RUN mkdir /build && cd /build \
 && /moab/configure \
      --enable-optimize \
      --disable-static \
      --enable-shared \
      --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial \
      --prefix=/usr/local \
 && make -j4 \
 && make install \
 && rm -Rf /build

# build and install DAGMC.  not sure why "make install" also copies 64MB of
# test scripts and data, let's clean those up
RUN mkdir /build && cd /build \
 && cmake /DAGMC -DMOAB_DIR=/usr/local -DCMAKE_INSTALL_PREFIX=/usr/local \
 && make -j4 all test \
 && make install \
 && rm -Rf /build /usr/local/tests

# build and install overlap_checker
RUN mkdir /build && cd /build \
 && cmake /overlap_checker -DCMAKE_INSTALL_PREFIX=/usr/local \
 && make -j4 \
 && make install \
 && rm -Rf /build

# build and install occ_faceter
RUN mkdir /build && cd /build \
 && cmake /occ_faceter -DCMAKE_BUILD_TYPE=Release -DMOAB_ROOT=/usr/local -DCMAKE_INSTALL_PREFIX=/usr/local \
 && make -j4 \
 && make install \
 && rm -Rf /build

# resulting image
FROM ubuntu_with_occ

# copy binaries that were built above
COPY --from=inner /usr/local/ /usr/local/

# also copy source over
COPY --from=inner /occ_faceter/ /occ_faceter/
COPY --from=inner /overlap_checker /overlap_checker/
