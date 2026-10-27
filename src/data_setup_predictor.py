import torch
from torch.utils.data import DataLoader, Dataset, random_split


class CustomDataset(Dataset):
    def __init__(self, data, target):
        self.data = torch.tensor(data, dtype=torch.float32)
        self.target = torch.tensor(target, dtype=torch.long)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        sample_x = self.data[index]
        sample_y = self.target[index]
        return sample_x, sample_y


DATA_MODE = "separate"
BATCH_SIZE = 32
NUM_WORKERS = 4


if DATA_MODE == "separate":
    X_rep = torch.load(
        "../data/embeddings/ori_predictor_separate/X_rep.pt"
    )
    y_rep = torch.load(
        "../data/embeddings/ori_predictor_separate/y_rep.pt"
    )

    X_relax = torch.load(
        "../data/embeddings/ori_predictor_separate/X_relax.pt"
    )
    y_relax = torch.load(
        "../data/embeddings/ori_predictor_separate/y_relax.pt"
    )

    X = torch.cat([X_rep, X_relax], dim=0)
    y = torch.cat([y_rep, y_relax], dim=0)

    full_dataset = CustomDataset(X, y)

    train_size = int(0.8 * len(full_dataset))
    val_size = len(full_dataset) - train_size

    print(f"Training Set: {train_size}")
    print(f"Validation Set: {val_size}")

    train_dataset, val_dataset = random_split(
        full_dataset,
        [train_size, val_size],
    )

    x_val_dataset = torch.stack(
        [val_dataset[i][0] for i in range(len(val_dataset))]
    )
    y_val_dataset = torch.stack(
        [val_dataset[i][1] for i in range(len(val_dataset))]
    )

    train_dataloader = DataLoader(
        dataset=train_dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=NUM_WORKERS,
        pin_memory=True,
    )

    print(x_val_dataset.shape)


elif DATA_MODE == "combined":
    X_train = torch.load(
        "../data/embeddings/ori_predictor_combined/combined_X_train.pt"
    )
    y_train = torch.load(
        "../data/embeddings/ori_predictor_combined/combined_y_train.pt"
    )
    X_val_rep = torch.load(
        "../data/embeddings/ori_predictor_combined/X_val_rep.pt"
    )
    y_val_rep = torch.load(
        "../data/embeddings/ori_predictor_combined/y_val_rep.pt"
    )
    X_val_relax = torch.load(
        "../data/embeddings/ori_predictor_combined/X_val_relax.pt"
    )
    y_val_relax = torch.load(
        "../data/embeddings/ori_predictor_combined/y_val_relax.pt"
    )

    combined_train_dataset = CustomDataset(X_train, y_train)

    train_dataloader = DataLoader(
        dataset=combined_train_dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=NUM_WORKERS,
        pin_memory=True,
    )

    x_val_dataset = torch.cat([X_val_rep, X_val_relax], dim=0)
    y_val_dataset = torch.cat([y_val_rep, y_val_relax], dim=0)

    perm_val = torch.randperm(len(x_val_dataset))
    x_val_dataset = x_val_dataset[perm_val]
    y_val_dataset = y_val_dataset[perm_val]

    print(X_train.shape)


else:
    raise ValueError(
        f"Unknown DATA_MODE: {DATA_MODE!r}. Use 'separate' or 'combined'."
    )