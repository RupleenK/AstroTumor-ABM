package OnLattice.base;

/**
 * The CellType class encapsulates the biological type and state of a cell.
 */
public class CellType {

    // Enum for the basic cell type.
    public enum celltype { TUMOR, ASTROCYTE, NEUTRAL }

    // Enum for the state of tumor cells.
    public enum TumorState { NORMAL }

    // Enum for the state of astrocyte cells.
    public enum AstrocyteState { ANTI_METASTATIC, PRO_METASTATIC }

    // The primary cell type.
    public celltype type;

    // The state of a tumor cell (if applicable).
    public TumorState tumorState;

    // The state of an astrocyte (if applicable).
    public AstrocyteState astrocyteState;

   // Constructs a new CellType with the given primary cell type.
    public CellType(celltype type) {
        this.type = type;
        // Initialize tumor or astrocyte state based on the cell type.
        switch (type) {
            case TUMOR:
                this.tumorState = TumorState.NORMAL; // Default tumor state.
                break;
            case ASTROCYTE:
                this.astrocyteState = AstrocyteState.ANTI_METASTATIC; // Default astrocyte state.
                break;
            default:
                // NEUTRAL cells do not require an additional state.
        }
    }

    // Returns the primary cell type.
    public celltype getCellTypeEnum() {
        return type;
    }

    // Returns the astrocyte state.
    public AstrocyteState getAstrocyteState() {
        return astrocyteState;
    }

    // Sets the primary cell type.
    public void setType(celltype newType) {
        this.type = newType;
    }

    // Sets the astrocyte state.
    public void setAstrocyteState(AstrocyteState newAstroState) {
        this.astrocyteState = newAstroState;
    }
}
