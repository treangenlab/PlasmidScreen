from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

from plasmidScreen.lib.funcs import get_default_db_path

custom_taxonomy = [
    (1, 1, "no rank", "root"),
    (2, 1, "superkingdom", "Cellular Organisms"),
    (99999, 2, "Species", "Plasmids"),
    (1012, 99999, "no rank", "Natural"),
    (1001, 99999, "no rank", "Synthetic"),
    (562, 2, "no rank", "Escherichia coli")
]
default_kmer_size = 21


@dataclass
class MetaMDBGConfig:
    """
    MetaMDBG configuration object for LR assembly
    """
    ont_mode: bool = True
    pacbio_mode: bool = False
    max_threads: int = 40
    kmer_size: int = 15


@dataclass
class MegahitConfig:
    """
    Megahit Configuration object for SR assembly
    """
    paired = True
    max_threads = 40


@dataclass
class KrakenConfig:
    """
    Configuration object.
    defined fields act as defaults.
    """
    # 1. Define your parameters and their default values here
    db_path: Path = get_default_db_path("KrakenEng", "default_storage.db")
    kmer_size: int = 21
    max_threads: int = 40
    masking: bool = True
    output_format: str = "csv"
    # must be equal to kmer_size for expected behavior to have kraken2 behave like kraken 1. As well as
    # minimizer-spaces=0
    minimizer_len: int = 21
    minimizer_spaces: int = 0
    load_factor = 0.7

    @classmethod
    def from_user_input(cls, user_overrides: dict[str, Any]) -> "AppConfig":
        """
        Creates a config object. Uses user_overrides where provided,
        falls back to class defaults otherwise.

        Safely ignores keys in user_overrides that don't match Config fields.
        """
        # Get the list of valid field names for this class
        valid_fields = cls.__annotations__.keys()

        # Keep only the keys that actually exist in our class and are not None
        # (Assuming None means "user didn't specify, use default")
        clean_overrides = {
            k: v for k, v in user_overrides.items()
            if k in valid_fields and v is not None
        }

        # Unpack the dictionary into the class constructor
        return cls(**clean_overrides)

    def to_dict(self) -> dict:
        """Returns the config as a dictionary (useful for logging/saving)."""
        return asdict(self)
