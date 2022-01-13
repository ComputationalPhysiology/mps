from PyInstaller.utils.hooks import collect_dynamic_libs

binaries = collect_dynamic_libs("ap_features")
if len(binaries) == 0:
    import ap_features

    binaries = [(ap_features.lib.lib._name, "ap_features")]
