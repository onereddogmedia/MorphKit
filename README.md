## MorphKit

Cut down version of [SpectMorph](https://github.com/swesterfeld/spectmorph)

### Compiling MorphKit:

Open MorphKit.xcodeproj with Xcode >11.3.1, build the MorphKit target Release configuration.

### Compiling FFTW

Use the included scripts to build the libraries, e.g.

```
cd fftw-3.3.9
./fftw-macos.sh
./fftw-ios.sh
```

Ensure you have appropriate Xcode command line tools installed. Adjust the variables as needed, e.g. XCODE_TOOLCHAINS

