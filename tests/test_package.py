import melscrna


def test_package_imports_and_has_version():
    assert melscrna.__version__ == "0.1.0"


def test_core_dependencies_import():
    import anndata  # noqa: F401
    import infercnvpy  # noqa: F401
    import scanpy  # noqa: F401
