from distutils.core import setup, Extension
extension_mod = Extension(
    "_crop",
    libraries=[
        "LightEnv",
    ],
    library_dirs=[
        "../../../../Lightenv"
    ],
    sources=[
        "crop.cc",
        "../crop/controller.cpp",
        "../crop/timer.cpp",
        "../crop/development.cpp",
        "../crop/plant.cpp",
        "../crop/initinfo.cpp",
        "../crop/radiation.cpp",
        "../crop/nodalunit.cpp",
        "../crop/gas_exchange.cpp",
        "../crop/organ.cpp",
        "../crop/ear.cpp",
        "../crop/thermaltime.cpp",
        "../crop/leaf.cpp",
        "../crop/stem.cpp",
        "../crop/roots.cpp",
        "../crop/cpmdate.cpp",
    ],
)
setup(name="crop", ext_modules=[extension_mod])
