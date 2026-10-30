"""Figure export shared by verification and results scripts."""
from pathlib import Path
import matplotlib.pyplot as plt

FIGURES_DIR = Path(__file__).resolve().parents[1] / 'figures'

def save_figure(name):
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / name
    plt.savefig(path, dpi=170, bbox_inches='tight')
    plt.close()
    print(f'Saved {path.relative_to(FIGURES_DIR.parent)}')
