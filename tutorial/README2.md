### 2.1. Installing demuxlet ###

#### 2.1.0. MAC NOTES ####

You will need to install some additional packages on a Mac.

First, you will need XCode CommandLineTools

Second, it's good to install homebrew:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Then, we need to install autoconf, automake, and libtool

```
$ brew install autoconf automake libtool
```

We will also need a new clang since Apple's clang doesn't support OpenMPI

```
$ brew install llvm
$ export CC=/usr/local/Cellar/llvm/4.0.0_1/bin/clang
$ export CXX=/usr/local/Cellar/llvm/4.0.0_1/bin/clang++ 
$ export LDFLAGS=-L/usr/local/Cellar/llvm/4.0.0_1/lib/
```

The CRAM format may use LZMA2 compression, which is implemented in HTSlib
by using compression routines from liblzma <http://tukaani.org/xz/>.

Building HTSlib requires liblzma development files to be installed on the
build machine; you may need to ensure a package such as liblzma-dev (on Debian
or Ubuntu Linux), xz-devel (on RPM-based Linux distributions or Cygwin), or
xz (via Homebrew on macOS) is installed; or build XZ Utils from source. Assuming we are on a Mac:

```
$ brew install xz
```

#### 2.1.1. Install htslib ####
```
$ git clone https://github.com/samtools/htslib.git
$ cd path/to/htslib
$ autoheader
$ autoconf
$ ./configure # optional, --prefix=/path/file
$ make
$ make install # optional, DESTDIR=/path/file
```

#### 2.1.2. Install demuxlet ####

demuxlet and htslib should be installed in same directory

```
$ git clone https://github.com/statgen/demuxlet.git
$ cd /path/to/demuxlet
$ mkdir m4
$ autoreconf -vfi
           
$ ./configure # optional, --prefix=/path/file
$ make
$ make install # optional, DESTDIR=/path/file
```
