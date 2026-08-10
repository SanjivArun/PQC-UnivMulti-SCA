#!/usr/bin/env python3
"""Generate a table image (PNG) comparing the Kyber NTT TVLA datasets."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── DATA ────────────────────────────────────────────────────────────────────
# Each inner list/tuple is one row, in the exact column order.
TABLE_COLUMNS = [
    "Dataset",
    "Implementation",
    "Traces (Total)",
    "Trace Length",
    "Purpose",
]

TABLE_DATA = [
    [
        "Unmasked",
        "ML-KEM NTT",
        "512",
        "93,142 samples",
        "TVLA",
    ],
    [
        "Masked",
        "ML-KEM NTT",
        "512",
        "168,352 samples",
        "TVLA",
    ],
]


# ─── RENDERING ──────────────────────────────────────────────────────────────
def save_table_png(
    path: str = "table.png",
    columns: list[str] | None = None,
    data: list[list] | None = None,
    nrows: int = 12,
    ncols: int = 5,
) -> None:
    columns = columns or TABLE_COLUMNS
    data = data or TABLE_DATA

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis("tight")
    ax.axis("off")

    table = ax.table(
        cellText=data,
        colLabels=columns,
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    table.scale(1.0, 1.9)

    # header
    for j in range(ncols):
        cell = table[0, j]
        cell.set_facecolor("#4472C4")
        cell.set_text_props(color="white", weight="bold", ha="center")

    # body rows
    for i in range(1, len(data) + 1):
        bg = "#D9E1F2" if i % 2 == 0 else "#FFFFFF"
        for j in range(ncols):
            cell = table[i, j]
            cell.set_facecolor(bg)
            cell.set_text_props(ha="center")

    plt.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ─── ENTRY POINT ────────────────────────────────────────────────────────────
def main() -> None:
    save_table_png()
    print("Saved to table.png")


if __name__ == "__main__":
    main()
