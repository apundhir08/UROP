import os, urllib.parse
def tiledb_cfg():
    p = urllib.parse.urlparse(os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy") or "")
    cfg = {}
    if p.hostname:
        cfg = {"vfs.s3.proxy_host": p.hostname,
               "vfs.s3.proxy_port": str(p.port),
               "vfs.s3.proxy_scheme": p.scheme or "http"}
    return cfg