# -*- mode: python ; coding: utf-8 -*-

a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=[],
    datas=[
        ('db/celestial_bodies.db', 'db'),
        ('kernels/de440.bsp', 'kernels'),
        ('kernels/naif0012.tls', 'kernels'),
        ('kernels/pck00010.tpc', 'kernels'),
	('workingFiles/astroScanLogo.png', 'workingFiles'),
    ],
    hiddenimports=['db', 'ephemerides', 'lambert', 'utils', 'plots'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='AstroScan',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['workingFiles\\icon.ico'],
)
