import torch
import torch.nn as nn

from data import train_dataloader, x_val_dataset, y_val_dataset
from losses import PositionWeightedCrossEntropy
from models import MultiHeadFourClassClassiferNN


SEED = 3242
EPOCHS = 30
LEARNING_RATE = 0.001
WEIGHT_DECAY = 1e-4

MODEL_CLASS = MultiHeadFourClassClassiferNN
MODEL_ARGS = (383, 64, 64)

#LOSS_FN = PositionWeightedCrossEntropy(rep_position_weights).to(device)
LOSS_FN = nn.CrossEntropyLoss


torch.manual_seed(SEED)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = MODEL_CLASS(*MODEL_ARGS).to(device)
loss_fn = LOSS_FN()
optimiser = torch.optim.Adam(
    model.parameters(),
    lr=LEARNING_RATE,
    weight_decay=WEIGHT_DECAY,
)


def accuracy_fn(y_true, y_pred):
    correct_pred = torch.eq(y_true, y_pred).sum().item()
    total_pred = y_true.numel()
    return (correct_pred / total_pred) * 100


acc = []
acc_val = []
loss_train = []
loss_val = []

for epoch in range(EPOCHS):
    model.train()

    running_loss = 0.0
    running_acc = 0.0

    for inputs, targets in train_dataloader:
        inputs = inputs.to(device, non_blocking=True)
        targets = targets.to(device, non_blocking=True)

        logits = torch.stack(model(inputs), dim=0).permute(1, 2, 0)
        predictions = torch.argmax(logits, dim=1)

        loss = loss_fn(logits, targets)

        optimiser.zero_grad()
        loss.backward()
        optimiser.step()

        running_loss += loss.item()
        running_acc += accuracy_fn(targets, predictions)

    epoch_loss = running_loss / len(train_dataloader)
    epoch_acc = running_acc / len(train_dataloader)

    model.eval()

    with torch.inference_mode():
        val_inputs = x_val_dataset.to(device, non_blocking=True)
        val_targets = y_val_dataset.to(device, non_blocking=True)

        val_logits = torch.stack(model(val_inputs), dim=0).permute(1, 2, 0)
        val_predictions = torch.argmax(val_logits, dim=1)

        val_loss = loss_fn(val_logits, val_targets)
        val_accuracy = accuracy_fn(val_targets, val_predictions)

    acc.append(epoch_acc)
    acc_val.append(val_accuracy)
    loss_train.append(epoch_loss)
    loss_val.append(val_loss.item())

    print(
        f"Epoch: {epoch + 1:02d}/{EPOCHS} | "
        f"Training Loss: {epoch_loss:.4f} | "
        f"Training Accuracy: {epoch_acc:.4f} | "
        f"Validation Loss: {val_loss.item():.4f} | "
        f"Validation Accuracy: {val_accuracy:.4f}"
    )