import itertools
import math

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec

# Note this data is already sorted
df = pd.read_csv("data/allCommonsLabeled.csv")


def matrixPlotter(df, Type, phenoClass):
    labelDictionary = pd.DataFrame(
        {
            "index": list(range(1, 15)),
            "name": [
                "DMN",
                "VIS",
                "FP",
                "DAN",
                "VAN",
                "SAL",
                "CO",
                "SMD",
                "SML",
                "AUD",
                "Tpole",
                "MTL",
                "PMN",
                "PON",
            ],
        }
    )
    if phenoClass == "Gordon":
        labels = pd.read_csv("brainTemplate/gordon_modules.csv")
    else:
        labels = pd.read_csv("brainTemplate/GRP1_template_parcel.csv")

    connectionOrder = (
        df.query("Type == 'twin'")
        .groupby(["phenoClass", "connection"])["C"]
        .mean()
        .sort_values(ascending=False)
        .reset_index()
    )

    connectionOrder = connectionOrder.loc[
        (connectionOrder.phenoClass == phenoClass),
    ].reset_index(drop=True)
    labels.columns = ["labels"]
    #    labels = labels.sort_values("labels")
    labels["name"] = [labelDictionary.name[idx - 1] for idx in labels.labels]
    axescounts = labels.labels.value_counts().sort_index()
    m = labels.shape[0]
    mat = np.zeros((m, m))
    # insert herit column into lower triangle of mat
    subdata = df[(df.phenoClass == phenoClass) & (df.Type == Type)]
    mat[np.tril_indices_from(mat, k=-1)] = subdata.C
    mat = mat + mat.T
    # mat[np.nan(mat)] = 0.00001
    # column wise mean without nan
    pd.DataFrame(np.nanmedian(mat, axis=1)).to_csv(
        f"formattedData/roiMedianHerit{phenoClass}_{Type}.csv"
    )

    # gordonNames = labels.name.unique()
    labels = labels.labels

    # Create labels for manhattan plot
    manLabels = labels.reset_index(drop=False)
    manLabels.columns = ["index", "name"]
    manLabelsMat = np.empty((m, m), object)

    for i in range(m):
        for j in range(m):
            manLabelsMat[i, j] = str(manLabels.name[i]) + "-" + str(manLabels.name[j])
    subdata.connection = pd.Categorical(
        subdata.connection, categories=connectionOrder.connection, ordered=True
    )
    subdata = (
        subdata.sort_values("connection").reset_index(drop=True).reset_index(drop=False)
    )

    # subdata[["first", "second"]] = subdata["connection"].str.split("-", expand=True)

    # # fill column first using the value to look up from labelDictionary
    # for row in range(len(subdata)):
    #     subdata["first"][row] = labelDictionary.iloc[int(subdata["first"][row]) - 1, 1]
    #     subdata["second"][row] = labelDictionary.iloc[
    #         int(subdata["second"][row]) - 1, 1
    #     ]

    # subdata["connection"] = subdata["first"] + "-" + subdata["second"]

    # reorder the manhattan groups by their sizes
    # sorted_groups = subdata.groupby("connection").size().sort_values(ascending=False)

    # make the "connection" in subdata a categorical and use sorted groups for order
    # subdata["connection"] = pd.Categorical(
    #     subdata["connection"], categories=sorted_groups.index
    # )
    # subdata = subdata.sort_values("connection").reset_index(drop=True)
    # subdata["withinNet"] = subdata["first"] == subdata["second"]
    # five number summary of h2 column grouped by withinNet column

    # print("############################################################")
    # print("############################################################")
    # print("############################################################")
    # print(Type, phenoClass)
    # print(subdata[["h2", "withinNet"]].groupby("withinNet").describe().T)
    # print("############################################################")
    # print("############################################################")
    # print("############################################################")

    manhattanDF_grouped = subdata.groupby("connection")

    # plot the grouped mean of the "signif" column in the manhattanDF_grouped dataframe
    # plot it as a barplot with the mean as y and the index as the x

    # label the top 20
    groupsizes = np.array(
        [group.shape[0] for num, (name, group) in enumerate(manhattanDF_grouped)]
    )
    # numbBigGroups = sum(groupsizes / sum(groupsizes) > 0.015)

    nlabels = 30
    pROIs = 0.005
    colors = np.fromiter(mcolors.TABLEAU_COLORS.values(), dtype="<U7")
    colors = np.tile(colors, math.ceil(nlabels / len(colors)))[range(nlabels)]
    colors = np.append(colors, np.repeat("grey", len(groupsizes) - nlabels))

    # fig, ax = plt.subplots(2, figsize=(3.5, 2.7))
    # fig, ax = plt.subplots(2)
    # fig.subplots_adjust(hspace=0)
    fig, ax = plt.subplots(1, figsize=(3.5, 1.15))
    # ax.set_title(f"{phenoClass} {Type} {m} ROIs")

    plotdata = manhattanDF_grouped["signif"].mean()
    colored = plotdata[0:nlabels]
    colored["others"] = 0
    # only select the first 30 grouped means
    # plotdata = plotdata[0:30]
    ax.bar(x=plotdata[0:nlabels].index, height=colored[0:nlabels], color=colors)
    plt.ylim([0, 0.35])
    # remove last ytick from ax
    ax.get_yticklabels()[-1].set_visible(False)
    ax.bar(
        x=np.linspace(
            start=nlabels, stop=nlabels + (nlabels / 10), num=len(plotdata[nlabels:])
        ),
        height=np.flip(plotdata[nlabels:]),
        color="grey",
        width=0.8 / 200,
    )
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        size=4.2,
        fontdict={"horizontalalignment": "right", "weight": "bold"},
    )
    ax.set_yticklabels(ax.get_yticklabels(), size=8)
    # ax.boxplot(plotdata[nlabels:], positions=[nlabels + ], tick_labels = ["other"])
    # ax.boxplot(plotdata[nlabels:], positions=[nlabels], tick_labels = ["other"])
    ax.set_ylabel("Prop heritable", size=8)

    plt.savefig(
        f"figures/fcs/manhattan{phenoClass}_{Type}_prop.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
    )
    minSigh2 = np.min(subdata[subdata.signif].C)

    fig, ax = plt.subplots(1, figsize=(3.5, 1.15))
    # manhattan plot
    x_labels = []
    x_labelsizes = []
    x_labels_pos = []
    divider = None
    for num, (name, group) in enumerate(manhattanDF_grouped):

        if divider is None:
            if num > (nlabels - 1):
                # if group.shape[0] / subdata.shape[0] <= pROIs:
                divider = min(group.index)
                print(f"divider: {divider}")
        if num > (nlabels - 1):
            # if group.shape[0] / subdata.shape[0] <= pROIs:
            group["index"] = (group["index"] - divider) / 100 + divider

        # for the following groups plot their indices as
        # (group.index - divider) / 10 + divider
        # make a vector of black and lightgrey based on whether
        # column signif is 0 or 1
        # colorVec = np.where(group.signif == 1, "black", "lightgrey")
        if num > (nlabels - 1):
            alphaVec = 0.1
        else:
            # alphaVec = np.where(group.signif == 1, 0.5, 0.1)
            alphaVec = 1
        group.plot(
            kind="scatter",
            x="index",
            y="Common",
            color=colors[num],
            # color="black",
            ax=ax,
            s=1,
            # alpha=alphaVec,
            zorder=10,
        )

        # dotted line to vizualize proportion
        if num < (nlabels):
            # if group.shape[0] / subdata.shape[0] > pROIs:
            # plot the group median as well
            # ax.plot(
            #     group["index"],
            #     # find the median where signif is 1
            #     np.repeat(
            #         # group["h2"].where(group["signif"] == 1).dropna().median(),
            #         # group["h2"].median(),
            #         group["signif"].mean(),
            #         len(group["index"]),
            #     ),
            #     # np.repeat(group["h2"].median(), len(group["index"])),
            #     zorder=10,
            #     color="black",
            #     linestyle="dashed",
            #     linewidth=1,
            # )
            # set background color locally
            # ax.axhspan(
            #     min(group.index),
            #     max(group.index) + 1,
            #     # facecolor="0.2",
            #     facecolor=colors[num % len(colors)],
            #     alpha=0.5,
            #     zorder=5,
            # )
            if group.shape[0] / subdata.shape[0] > pROIs:
                x_labels.append(name)
                x_labelsizes.append(group.shape[0] / subdata.shape[0])
                x_labels_pos.append(
                    (
                        group["index"].iloc[-1]
                        - (group["index"].iloc[-1] - group["index"].iloc[0]) / 2
                    )
                )
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(
        x_labels,
        rotation=45,
        size=4.2,
        fontdict={"horizontalalignment": "right", "weight": "bold"},
    )
    ax.set_xlabel("")
    ax.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        # bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        # labelbottom=False,
    )  # labels along the bottom edge are off
    # plot minimum signif h2 as a dotted line
    # ax.plot(
    #     subdata["index"],
    #     # find the median where signif is 1
    #     np.repeat(
    #         # group["h2"].where(group["signif"] == 1).dropna().median(),
    #         # group["h2"].median(),
    #         minSigh2,
    #         len(subdata["index"]),
    #     ),
    #     zorder=10,
    #     color="black",
    #     linestyle="dashed",
    #     linewidth=1,
    # )

    # set axis limits
    ax.set_xlim([0, max(group["index"])])
    ax.set_ylim([0, 1])

    ax.set_ylabel(r"Heritability ($h^2$)", size=8)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1], size=8)
    # plot vertical line for the mean over all h2 values
    # ax.axvline(x=subdata["h2"].mean(), color="red", linestyle="dashed", linewidth=1)
    # adjust this so that the figure is 2.5 in x 1.3 in
    plt.savefig(
        f"figures/fcs/manhattan{phenoClass}_{Type}_mag.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
    )
    # plt.savefig(
    #     f"figures/fcs/manhattan{phenoClass}_{Type}.png",
    #     dpi=300,
    #     bbox_inches="tight",
    #     pad_inches=0,
    # )

    # gs_kw = dict(
    #     width_ratios=axescounts.values,
    #     height_ratios=axescounts.values,
    #     hspace=0.1,
    #     wspace=0.1,
    # )
    # fig, axes = plt.subplots(
    #     len(np.unique(labels)), len(np.unique(labels)), gridspec_kw=gs_kw
    # )

    # for r, row in enumerate(axes):
    #     for c, ax in enumerate(row):
    #         # replace nas with zer in mat
    #         mat = np.where(mat == np.nan, 0.0001, mat)
    #         ax.imshow(
    #             mat[labels == r + 1, :][:, labels == c + 1],
    #             vmin=0,
    #             vmax=max(subdata.h2),
    #             cmap="hot",
    #             aspect="auto",
    #         )
    #         ax.set(xlabel="", ylabel="")
    #         ax.set_xticklabels([])
    #         ax.set_xticks([])
    #         ax.set_yticklabels([])
    #         ax.set_yticks([])
    #         if r == 0:
    #             ax.set_xlabel(labelDictionary.name[c], rotation=45)
    #             ax.xaxis.set_label_position("top")
    #         if c == 0:
    #             ax.set_ylabel(labelDictionary.name[r], rotation=45)

    # # add a colormap legend outside of the subplots
    # # Create ScalarMappable object for colorbar
    # sm = plt.cm.ScalarMappable(
    #     cmap="hot", norm=plt.Normalize(vmin=0, vmax=max(subdata.h2))
    # )
    # sm.set_array([])  # dummy array for the colorbar

    # # Add colorbar
    # plt.colorbar(sm, ax=axes.ravel().tolist())

    # fig.savefig(f"figures/fcs/freeHeatmap{phenoClass}_{Type}.png", dpi=300)

    # # save subdata to csv
    # subdata.to_csv(f"formattedData/subdata{phenoClass}_{Type}.csv")
    return ax


# figs = []
# matrixPlotter(df, Type="Twin", phenoClass="PFN")
for Type, phenoClass in itertools.product(["twin"], np.unique(df.phenoClass)):
    matrixPlotter(df, Type=Type, phenoClass=phenoClass)
