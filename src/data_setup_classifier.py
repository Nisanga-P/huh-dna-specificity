from pathlib import Path

import torch
from torch.utils.data import DataLoader, Dataset, random_split


class CustomDataset(Dataset):
    def __init__(self, data, target):
        self.data = torch.as_tensor(data, dtype=torch.float32)
        self.target = torch.as_tensor(target, dtype=torch.float32)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        sample_x = self.data[index]
        sample_y = self.target[index]
        return sample_x, sample_y


tensor_dir = Path("../data/embeddings/classifier_combined/")

X_train = torch.load(tensor_dir / "X_train.pt", weights_only=True)
y_train = torch.load(tensor_dir / "y_train.pt", weights_only=True)
X_val_rep = torch.load(tensor_dir / "X_val_rep.pt", weights_only=True)
y_val_rep = torch.load(tensor_dir / "y_val_rep.pt", weights_only=True)
X_val_relax = torch.load(tensor_dir / "X_val_relax.pt", weights_only=True)
y_val_relax = torch.load(tensor_dir / "y_val_relax.pt", weights_only=True)

print(X_train.shape)

full_dataset = CustomDataset(X_train, y_train)

train_size = int(0.8 * len(full_dataset))
test_size = len(full_dataset) - train_size

train_dataset, test_dataset = random_split(
    full_dataset,
    [train_size, test_size],
)

count_p = sum(1 for _, target in train_dataset if target.item() == 1)
print(f"(Train) Positive: {count_p}")
print(f"(Train) Negative: {len(train_dataset) - count_p}")

count = sum(1 for _, target in test_dataset if target.item() == 1)
print(f"(Test) Positive: {count}")
print(f"(Test) Negative: {len(test_dataset) - count}")

x_test_dataset = torch.stack(
    [test_dataset[i][0] for i in range(len(test_dataset))]
)
y_test_dataset = torch.stack(
    [test_dataset[i][1] for i in range(len(test_dataset))]
)

train_dataloader = DataLoader(
    dataset=train_dataset,
    batch_size=32,
    shuffle=True,
    num_workers=4,
    pin_memory=True,
)