import json
import os
import sys
from pathlib import Path


def parse_backbone_map(file_path):
    """
    Parses a JSON file to create a dictionary mapping 'id' to 'backbone'.

    Args:
        file_path (str): The path to the JSON file.

    Returns:
        dict: A dictionary {id_value: backbone_value}.
    """
    result_map = {}

    try:
        with open(file_path, 'r') as f:
            data = json.load(f)

        plasmids_list = data["plasmids"]
        for plasmid in plasmids_list:
            result_map[plasmid["id"]] = plasmid["cloning"]["backbone"]

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except json.JSONDecodeError:
        print(f"Error: Failed to decode JSON from '{file_path}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return result_map


def get_default_db_path(app_name: str, db_filename: str = "app.db") -> str:
    """
    Determines the platform-specific default data directory and returns
    the full path to the database file. Ensures the directory exists.

    Args:
        app_name: The name of your application (used for the folder name).
        db_filename: The name of the database file.

    Returns:
        Path: The full path to the database file.
    """
    home = Path.home()

    if sys.platform == "win32":
        # Windows: %LOCALAPPDATA% -> C:\Users\Name\AppData\Local\AppName
        # Fallback to %APPDATA% if LOCALAPPDATA is missing
        base_dir = Path(os.environ.get("LOCALAPPDATA", home / "AppData" / "Local"))
        data_dir = base_dir / app_name

    elif sys.platform == "darwin":
        # macOS: ~/Library/Application Support/AppName
        data_dir = home / "Library" / "Application Support" / app_name

    else:
        # Linux/Unix: ~/.local/share/app_name (XDG_DATA_HOME default)
        # Check XDG_DATA_HOME environment variable first
        xdg_data = os.environ.get("XDG_DATA_HOME")
        if xdg_data:
            data_dir = Path(xdg_data) / app_name
        else:
            data_dir = home / ".local" / "share" / app_name

    # Create the directory if it doesn't exist
    data_dir.mkdir(parents=True, exist_ok=True)

    return str(data_dir)# / db_filename

