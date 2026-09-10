#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sb


def load_batch_data(paths):
    """Load per-batch CSVs produced by report_stats.py --batch.

    Returns a dict mapping descriptor -> DataFrame(index=stats, cols=batch_values),
    ordered by file creation time.  Descriptor is inferred from the filename
    report_stats_{descriptor}_batch.csv.
    """
    batch_data = {}
    for path in sorted(paths, key=os.path.getctime):
        name = os.path.basename(path)
        descriptor = name.removeprefix('report_stats_').removesuffix('_batch.csv')
        batch_data[descriptor] = pd.read_csv(path, header=0, index_col=0)
    return batch_data


def format_label(col):
    return col[0].upper() + col[1:].replace('_', ' ').replace('mread', 'read').replace('tscp', 'transcript')


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('csv_list', nargs='+')
    p.add_argument('--batch_csvs', nargs='+', default=None)
    return p


def main():
    args = build_parser().parse_args()

    csv_paths = sorted(args.csv_list, key=os.path.getctime)
    df = pd.concat([pd.read_csv(f, header=0, index_col=0) for f in csv_paths], axis=1).T

    batch_data = load_batch_data(args.batch_csvs) if args.batch_csvs else None

    df_total = df[[col for col in df.columns if col.startswith('total')]].sort_values(by='total_nuclei', axis=0)

    # ------------------------------------------------------------------ bar chart
    fig, axs = plt.subplots(len(df_total.columns), 1, figsize=(8, 3 * len(df_total.columns)))

    if batch_data:
        batch_vals = list(next(iter(batch_data.values())).columns)
        n_groups = 1 + len(batch_vals)
        bar_h = 0.7 / n_groups
        palette = ['steelblue'] + list(sb.color_palette('Set2', n_colors=len(batch_vals)))
        steps = list(df_total.index)
        y_base = np.arange(len(steps))

        def draw_grouped_bars(ax, col, is_nuclei=False):
            top_y = y_base + (n_groups - 1) / 2 * bar_h
            # overall
            vals = df_total[col].values.astype(float)
            widths = vals if is_nuclei else np.log1p(vals)
            labels = [f'{v:,.0f}' if is_nuclei else f'{v:.2e}' for v in vals]
            p = ax.barh(top_y, widths, height=bar_h, color=palette[0], label='all')
            ax.bar_label(p, labels=labels, label_type='edge', padding=5)
            # one bar per batch
            for b_idx, bv in enumerate(batch_vals):
                raw = np.array([batch_data[s].loc[col, bv] if s in batch_data else np.nan
                                for s in steps], dtype=float)
                widths_b = np.nan_to_num(raw if is_nuclei else np.log1p(raw), nan=0)
                labels_b = [f'{v:,.0f}' if is_nuclei else f'{v:.2e}' if pd.notna(v) else ''
                            for v in raw]
                y = top_y - (b_idx + 1) * bar_h
                p = ax.barh(y, widths_b, height=bar_h, color=palette[b_idx + 1], label=bv)
                ax.bar_label(p, labels=labels_b, label_type='edge', padding=5)
            ax.set_yticks(y_base)
            ax.set_yticklabels(steps)
            ax.set_xlim(xmin=0)
            ax.spines[['top', 'right']].set_visible(False)

        for i, col in enumerate([x for x in df_total.columns if x != 'total_nuclei']):
            draw_grouped_bars(axs[i], col)
            axs[i].set_title(f'log1p {col.replace("_", " ")}', loc='left')
        draw_grouped_bars(axs[-1], 'total_nuclei', is_nuclei=True)
        axs[-1].set_title('total cells', loc='left')
        axs[0].legend(loc='lower right', fontsize=8)

    else:
        for i, col in enumerate([x for x in df_total.columns if x != 'total_nuclei']):
            data = df_total[col]
            for index, row in data.items():
                p = axs[i].barh(y=index, width=np.log1p(row))
                axs[i].bar_label(p, labels=[f'{row:.2e}'], label_type='edge', padding=5)
            axs[i].set_xlim(xmin=0)
            axs[i].set_title(f'log1p {col.replace("_", " ")}', loc='left')
            axs[i].spines[['top', 'right']].set_visible(False)

        for index, row in df_total['total_nuclei'].items():
            p = axs[-1].barh(y=index, width=row)
            axs[-1].bar_label(p, label_type='edge', fmt='{:,.0f}', padding=5)
        axs[-1].set_xlim(xmin=0)
        axs[-1].set_title('total cells', loc='left')
        axs[-1].spines[['top', 'right']].set_visible(False)

    fig.suptitle('Total counts after each filtering step', fontsize=20)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig("plt_filter_bars_summary.png", bbox_inches='tight', transparent=True, dpi=300)
    plt.close(fig)

    # ------------------------------------------------------------------ line chart
    df_median = df[[col for col in df.columns if col.startswith('median')]]

    fig, ax = plt.subplots()
    metric_colors = ['red', 'purple', 'blue']
    metric_linestyles = ['-', '--', ':']

    if batch_data:
        batch_vals = list(next(iter(batch_data.values())).columns)
        batch_palette = sb.color_palette('Set2', n_colors=len(batch_vals))
        steps = list(df_median.index)

        for i, col in enumerate(df_median.columns):
            lbl = format_label(col)
            ls = metric_linestyles[i % len(metric_linestyles)]
            # overall: metric colour, thicker line
            ax.plot(steps, df_median[col], ls=ls, lw=2, marker='o',
                    color=metric_colors[i], label=f'{lbl} (all)')
            # one line per batch: batch colour, same linestyle encodes metric
            for b_idx, bv in enumerate(batch_vals):
                y_vals = [batch_data[s].loc[col, bv] if s in batch_data else np.nan
                          for s in steps]
                ax.plot(steps, y_vals, ls=ls, lw=1, marker='.',
                        color=batch_palette[b_idx], label=f'{lbl} ({bv})')

        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8)

    else:
        for i, col in enumerate(df_median.columns):
            lbl = format_label(col)
            ax.plot(df_median.index, df_median[col], ls='--', lw=1, marker='o',
                    color=metric_colors[i], label=lbl)
            for x, y in zip(df_median.index, df_median[col]):
                ax.annotate(f'{y:,.0f}', xy=(x, y), textcoords='offset points',
                            xytext=(0, 6), ha='center', fontsize=10)
            ax.annotate(lbl, xy=(x, y), textcoords='offset points',
                        xytext=(6, -6), ha='left', fontsize=10, color=metric_colors[i])

    ax.set_ylim(ymin=0)
    ax.spines[['top', 'right']].set_visible(False)
    plt.suptitle('Median counts at each filtering step')
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig("plt_filter_lines_summary.png", bbox_inches='tight', transparent=True, dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
