#!/usr/bin/env python3

import os.path
import argparse
import traceback

import pandas as pd
import plotext


"""
Quickly plot data on terminal and Matplotlib interactive plot

-----------------------------------------
use "python3 plot_data.py -h" for usage
-----------------------------------------
NOTE: Uses WHITESPACES as delimiter by default, can be changed with --sep option

EXAMPLE: ./plot_data.py [--sep=','] [-x='<X_COL_NAME>'] [-y='<Y_COL_NAME>'] [-xl='<X_LABEL>'] [-yl='<Y_LABEL>'] [-t='PLOT TITLE'] [-o='<out-fig>'] data.csv
"""

COMMENT_TOKEN = "#"


def is_empty(s: str) -> bool:
    return s is None or s.strip() == ""


def print_full_line(char='-'):
    """
    Prints a line that fills the entire width of the terminal.
    """
    try:
        terminal_width = os.get_terminal_size().columns
        print(char * terminal_width)
    except OSError:
        # Handle cases where terminal size cannot be determined (e.g., not in a TTY)
        # print("Unable to determine terminal size. Printing a default line.")
        print(char * 80)  # Print a line of 80 characters as a fallback


def plot_data(data_file: str,
              delimiter_char: str = None,
              x_col_name: str = None,
              y_col_name: str = None,
              x_col_label: str = None,
              y_col_label: str = None,
              title: str = None,
              show_interactive_plot: bool = True,
              plot_terminal: bool = True,
              out_fig_file: str = None):

    sep = r"\s+"  # Defaults to whitespaces
    if not is_empty(delimiter_char):
        if len(delimiter_char) == 1:
            sep = delimiter_char
        else:
            print(f"ERROR: Invalid delimiter char \"{delimiter_char}\", Must be a single character !!")
            return

    if is_empty(title):
        title = os.path.splitext(data_file)[0]

    df = pd.read_csv(data_file, sep=sep, comment=COMMENT_TOKEN)
    if len(df.columns) == 0:
        print(f"WARNING: NO DATA TO PLOT. Input data file: \"{data_file}\"")
        return

    if len(df.columns) == 1:
        print(f"WARNING: Only 1 column in input data file: \"{data_file}\". Ignoring input column names...")
        x_col_name = "Index"
        x_col_label = "Index"
        y_col_name = df.columns[0]

        x_col = df.index
        y_col = df[y_col_name]
    else:
        if is_empty(x_col_name):
            x_col_name = df.columns[0]
            x_col_label = None

        if is_empty(y_col_name):
            y_col_name = df.columns[1]
            y_col_label = None

        x_col = df[x_col_name]
        y_col = df[y_col_name]

    # Plot on terminal
    if plot_terminal:
        plotext.theme("pro")

        plotext.plot(x_col, y_col)
        plotext.xlabel(x_col_name if is_empty(x_col_label) else x_col_label)
        plotext.ylabel(y_col_name if is_empty(y_col_label) else y_col_label)

        plotext.title(title)

        print("")
        print_full_line("#")
        print("")
        plotext.show()
        print("")
        print_full_line("#")
        print("")
    else:
        print("LOG: skipping TERMINAL PLOT")

    # Matplotlib Plot
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("ERROR: COuld not import Matplotlib!!")
        traceback.print_exc()
        return

    plt.plot(x_col, y_col)
    plt.xlabel(x_col_name if is_empty(x_col_label) else x_col_label)
    plt.ylabel(y_col_name if is_empty(y_col_label) else y_col_label)

    if is_empty(title):
        title = os.path.splitext(data_file)[0]
    plt.title(title)
    plt.ticklabel_format(axis='both', style='sci')

    # Save Figure
    if is_empty(out_fig_file):
        out_fig_file = os.path.splitext(data_file)[0] + '.pdf'
    plt.savefig(out_fig_file)

    # Interactive Plot
    if show_interactive_plot:
        plt.show()
    else:
        print("LOG: skipping INTERACTIVE PLOT")


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Quickly plot data on Terminal and Matplotlib.")

    argparser.add_argument("data_file", type=str, help="Input data file")
    argparser.add_argument("-s", "--sep", type=str, help=f"Separator character (defaults to whitespaces)")
    argparser.add_argument("-t", "--title", type=str, help="Plot title")
    argparser.add_argument("-x", "--x-col", type=str, help="X column to read from data file")
    argparser.add_argument("-y", "--y-col", type=str, help="Y column to read from data file")
    argparser.add_argument("-xl", "--x-label", type=str, help="X Label")
    argparser.add_argument("-yl", "--y-label", type=str, help="Y Label")
    argparser.add_argument("-o", "--out-fig", type=str,
                           help="Output figure file name with extension (.svg, .pdf, .png etc)")
    argparser.add_argument("-ni", "--no-interactive-plot", action="store_true", help="DO NOT show interactive plot")
    argparser.add_argument("-nt", "--no-terminal", action="store_true", help="DO NOT plot on terminal")

    args = argparser.parse_args()

    plot_data(data_file=args.data_file,
              delimiter_char=args.sep,
              x_col_name=args.x_col,
              y_col_name=args.y_col,
              x_col_label=args.x_label,
              y_col_label=args.y_label,
              title=args.title,
              show_interactive_plot=not args.no_interactive_plot,
              plot_terminal=not args.no_terminal,
              out_fig_file=args.out_fig)
