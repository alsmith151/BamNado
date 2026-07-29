"""Console entry point forwarding to the Rust clap CLI."""

import contextlib
import signal
import sys
import threading

from ._bamnado import _cli_main


@contextlib.contextmanager
def _default_sigint():
    """Restore OS-default SIGINT while the GIL is released by the Rust run."""
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    try:
        yield
    finally:
        signal.signal(signal.SIGINT, previous)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI and return its exit code.

    SIGINT uses the OS-default terminating behaviour during execution, matching
    the standalone binary rather than raising ``KeyboardInterrupt``.
    """
    args = sys.argv[1:] if argv is None else list(argv)
    with _default_sigint():
        return _cli_main(["bamnado", *args])
