
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict

@dataclass
class PipelineState:
    outputs: Dict[str, Path] = field(default_factory=dict)

    def register(self, key: str, path: Path):
        self.outputs[key] = Path(path)

    def get(self, key: str) -> Path:
        return self.outputs[key]
    
    def set(self, key, value):
        self.outputs[key] = value

