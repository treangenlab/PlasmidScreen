from __future__ import annotations

import tempfile

import os
import sys
from pathlib import Path


def get_default_db_path(app_name: str) -> str:
    """
    Determines the platform-specific default data directory and returns
    the full path to the database file. Ensures the directory exists.

    Args:
        app_name: The name of your application (used for the folder name).

    Returns:
        The full path to the data directory.
    """
    home = Path.home()

    if sys.platform == "win32":
        base_dir = Path(os.environ.get("LOCALAPPDATA", home / "AppData" / "Local"))
        data_dir = base_dir / app_name

    elif sys.platform == "darwin":
        data_dir = home / "Library" / "Application Support" / app_name

    else:
        xdg_data = os.environ.get("XDG_DATA_HOME")
        if xdg_data:
            data_dir = Path(xdg_data) / app_name
        else:
            data_dir = home / ".local" / "share" / app_name

    try:
        data_dir.mkdir(parents=True, exist_ok=True)

        # Test write permissions explicitly (handles cases where mkdir succeeds but write fails)
        if not os.access(data_dir, os.W_OK):
            raise PermissionError

    except (PermissionError, OSError):
        # 4. Fallback Strategy: Use the OS temporary directory
        # tempfile.gettempdir() automatically resolves to /tmp on Linux/macOS
        data_dir = Path(tempfile.gettempdir()) / app_name
        data_dir.mkdir(parents=True, exist_ok=True)

        # Alert the user/logs that data persistence is volatile
        print(
            f"Warning: Primary data directory unwritable. "
            f"Falling back to temporary storage: {data_dir}",
            file=sys.stderr
        )

    return str(data_dir)